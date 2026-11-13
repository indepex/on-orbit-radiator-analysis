import warnings
import sys
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun, GeocentricTrueEcliptic, SkyCoord
from astropy import units
from sgp4.earth_gravity import wgs84
from sgp4.api import Satrec
from datetime import timedelta
from scipy.optimize import minimize_scalar
import quaternion
from astropy.time import TimeDelta

from lib_utils import normalize


# -----------------------------------------------------------------------------------------
def earth_couple(dirs, dir_emission, toSat, SatX, SatY, SatZ, altitude):
    """
    Given satellite position vectors toSat, body axes SatX/Y/Z, and directional
    emission dir_emission for dirs, compute the emission intercepted by Earth.
    dirs         : (M, 3) directions in the satellite body frame
    dir_emission : (M,) emission associated with each direction
    toSat        : (N, 3) geocenter-to-satellite direction vectors
    SatX/Y/Z     : (N, 3) body-axis directions in the inertial frame

    return:
        couple : (N,) emission intercepted by Earth (space is also returned)
    """
    earth = normalize(-toSat)   # direction toward Earth from the satellite, (N, 3)
    e_body = np.column_stack([np.sum(SatX * earth, axis=1), np.sum(SatY * earth, axis=1),
                             np.sum(SatZ * earth, axis=1)])  # (N,3)
    Re = wgs84.radiusearthkm  # Earth radius [km]
    alpha = np.arcsin(Re / (Re + altitude))
    cos_alpha = np.cos(alpha)

    # cosang[n, m] = dot product of e_body[n] and dirs[m]
    cosang = e_body @ dirs.T     # (N,M)

    mask = cosang >= cos_alpha   # (N,M)

    space = np.sum((~mask) * dir_emission[None, :], axis=1)
    couple = np.sum(mask * dir_emission[None, :], axis=1)

    return couple, space


# -----------------------------------------------------------------------------------------
def load_outward(datafile):
    data = np.load(datafile)
    modelfile    = data['modelfile'].item()
    num_total    = data['num_rays_total'].item()
    total        = data['total_emission'].item()
    absorbed     = data['absorbed_emission'].item()
    Dirs         = data['Dirs']
    dir_emission = data['dir_emission']
    npoint       = data['Npoint']
    area         = data['Area']

    return modelfile, num_total, total, absorbed, Dirs, dir_emission, npoint, area


# -----------------------------------------------------------------------------------------
def get_attitude(orb, obj, npos, delay=0):
    """
    Compute satellite attitude over one orbital revolution.

    Returns
    -------
    toSat : ndarray (N, 3)
        Unit direction vectors from the geocenter to the satellite.
    theta : ndarray (N,)
        Orbital phase [deg], with zero at the start of observation.
    SatX, SatY, SatZ : ndarray (N, 3)
        Unit vectors along the satellite body axes.
    jdc : astropy.time.Time
        Central time of the observing half-orbit.
    toSun : ndarray (N, 3)
        Geocentric unit direction vectors toward the Sun.
    ranmid : float
        Roll correction at the midpoint of the observing interval [deg].
    time_jd : ndarray (N,)
        Julian dates of the samples (UTC).
    position : ndarray (N, 3)
        Satellite positions (TEME, km).
    velocity : ndarray (N, 3)
        Satellite velocities (TEME, km/s).

    Notes
    -----
    N = npos + 1.

    Satellite positions and velocities are in the SGP4 TEME frame.
    Sun directions are geocentric Astropy directions; target directions are in ICRS.
    The attitude calculation approximates these as a common inertial frame.
    No rigorous transformation between these frames is performed.
    """
    # Shift the observation start according to delay (typically 0, 1, or 2)
    jd = get_jd_from_sun_lon(obj['year'], obj['sun_longitude'] * units.deg)
    # Center of the observing half-orbit in the revolution containing that Julian date
    jdc = center_jd(jd, orb['satellite'], orb['Torb'], obj['toObj'])
    # Define phase zero at observation start, even when shifted by delay
    # theta = np.linspace(360.0 / npos * delay, 360.0 + 360.0 / npos * delay, npos + 1)
    theta = np.linspace(0.0, 360.0, npos + 1)   # npos + 1 samples, including the start of the next revolution

    # Adjust observation start time for delay
    # Generate time samples for one revolution
    with warnings.catch_warnings():  # Limit warning suppression to time calculations and SGP4
        warnings.simplefilter("ignore", ErfaWarning)
        warnings.simplefilter("ignore", AstropyWarning)
        dt = delay * orb["Torb"] / npos  # Observation-time offset due to delay
        obs_start_time = - orb['Torb'] / 4 + dt    # Start offset from a quarter orbit, adjusted for delay
        obs_next_time = orb['Torb'] * 3 / 4 + dt   # Start of the next revolution
        # jdlist = np.linspace(-orb['Torb'] / 4 + dt, 3 * orb['Torb'] / 4 + dt, npos + 1)
        jdlist = np.linspace(obs_start_time, obs_next_time, npos + 1)
        jds = jdc + TimeDelta(jdlist, format="sec")  # Observation times
        # Sun-direction array toSun
        toSun = normalize(get_sun(jds).cartesian.xyz.value.T)
        # Satellite positions
        e, Sat, Vel = orb["satellite"].sgp4_array(jds.jd1, jds.jd2)
        if np.any(e != 0):
            raise RuntimeError(f"SGP4 error codes: {np.unique(e)}")

    toSat = normalize(Sat)
    # Determine orbital rotation sense from LTAN
    if orb["ltan"] < 12:
        sign = -1
    else:
        sign = 1
    # Attitude during observation
    Zo, Yo, Xo, ranmid = get_obs_attitude(toSun[:npos // 2 + 1 - delay * 2], obj)
    # Attitude at the end of observation
    X0, Y0, Z0 = Xo[npos // 2 - delay * 2], Yo[npos // 2 - delay * 2], Zo[npos // 2 - delay * 2]
    # Attitude at the start of the next revolution
    Z1r, Y1r, X1r, dummy = get_obs_attitude(toSun[npos:npos + 1], obj)
    Z1, Y1, X1 = Z1r[0], Y1r[0], X1r[0]
    # Attitude outside the observing interval
    Xi, Yi, Zi = get_idle_attitude(X0, Y0, Z0, X1, Y1, Z1, toSun, toSat, npos, sign, delay)
    SatX = np.vstack([Xo, Xi])
    SatY = np.vstack([Yo, Yi])
    SatZ = np.vstack([Zo, Zi])

    return toSat, theta, SatX, SatY, SatZ, jdc, toSun, ranmid, jds.jd, Sat, Vel


# -----------------------------------------------------------------------------------------
def slerp_controlled(q0, q1, u, shortest=True):
    """
    shortest=True  : interpolate along the shorter arc
    shortest=False : interpolate along the longer arc
    """
    q0 = q0.normalized()
    q1 = q1.normalized()

    dot = q0.w * q1.w + q0.x * q1.x + q0.y * q1.y + q0.z * q1.z

    if shortest:
        if dot < 0.0:
            q1 = -q1
            dot = -dot
    else:
        if dot > 0.0:
            q1 = -q1
            dot = -dot

    dot = np.clip(dot, -1.0, 1.0)

    if dot > 0.9995:
        q = (1 - u) * q0 + u * q1
        return q.normalized()

    theta = np.arccos(dot)
    sin_t = np.sin(theta)

    a = np.sin((1 - u) * theta) / sin_t
    b = np.sin(u * theta) / sin_t
    q = a * q0 + b * q1
    return q.normalized()


# -----------------------------------------------------------------------------------------
def make_q01_and_interp_frame(X0, Y0, Z0, X1, Y1, Z1, use_shortest=True):
    """
    Return the quaternion rotating frame (X0, Y0, Z0) to (X1, Y1, Z1),
    and frame_at(u), which returns interpolated axes (Xt, Yt, Zt) for u in [0, 1].
    """
    A = np.column_stack([X0, Y0, Z0])
    B = np.column_stack([X1, Y1, Z1])
    # Rotation that maps axes of A to axes of B:  R = B * A^T
    R01 = B @ A.T
    q01 = quaternion.from_rotation_matrix(R01)

    # identity quaternion for interpolation start
    qI = np.quaternion(1.0, 0.0, 0.0, 0.0)

    def frame_at(u):
        q = slerp_controlled(qI, q01, u, shortest=use_shortest)
        Rt = quaternion.as_rotation_matrix(q)      # 3x3 rotation matrix
        Ft = Rt @ A                     # rotate t0-frame axes
        Xt, Yt, Zt = Ft[:, 0], Ft[:, 1], Ft[:, 2]
        return Xt, Yt, Zt

    return q01, frame_at


# -----------------------------------------------------------------------------------------
def get_idle_attitude(X0, Y0, Z0, X1, Y1, Z1, toSun, toSat, npos, sign, delay):
    """
    Attitude during the non-observing interval.
    Index npos//2 - 2 * delay is the last observation sample and is excluded.
    Its attitude is (X0, Y0, Z0).
    Index npos marks the next observation start with axes (X1, Y1, Z1); include it.
    """
    trans_angle = 10   # Use a 10-degree orbital phase interval for the transition.
    di = npos * trans_angle // 360  # Number of samples in the transition; chosen for integer division.
    """
    Over the non-observing interval, the satellite rotates by approximately 360 degrees.
    Except during transitions, the rotation axis and satellite X axis point toward the Sun.
    At the interval midpoint, the satellite, Sun, and geocenter define a plane; X points sunward.
    """
    idc = npos // 4 * 3 - 2 * delay  # Index of toSun at the midpoint of the non-observing interval
    Xic = toSun[idc]   # X-axis direction at the non-observing midpoint
    # Z points away from Earth within the satellite-Sun-geocenter plane
    Zic = normalize(toSat[idc] - np.dot(toSat[idc], toSun[idc]) * toSun[idc])
    Yic = normalize(np.cross(Zic, Xic))
    # Allocate attitude arrays for the non-observing interval
    X = np.empty((npos // 2 + 2 * delay + 1, 3))  # Last entry is the start of the next observation
    Y = np.empty((npos // 2 + 2 * delay + 1, 3))
    Z = np.empty((npos // 2 + 2 * delay + 1, 3))
    """
    Outside the initial and final transitions, rotate steadily about the sunward satellite X axis.
    The full non-observing interval spans 360 degrees over npos//2 + 2 * delay steps.
    """
    idcc = npos // 4 + delay
    for i in range(di, npos // 2 + 2 * delay  - di + 1):
        phi = sign * np.radians((i - idcc) * 360 / (npos / 2 + 2 * delay))
        qrot = quaternion.from_rotation_vector(Xic * phi)    # Rotation about X = toSun
        X[i] = Xic
        Y[i] = normalize(quaternion.as_rotation_matrix(qrot).dot(Yic))
        Z[i] = normalize(quaternion.as_rotation_matrix(qrot).dot(Zic))
        # print(f"idle {i} {i-idcc} phi {np.degrees(phi)}")

    """
    The initial transition starts from the final observing attitude (X0, Y0, Z0)
    and ends at the first steady-rotation attitude (X[di], Y[di], Z[di]).
    """
    # Quaternion and frame interpolation function for the transition
    q1a, frame1 = make_q01_and_interp_frame(X0, Y0, Z0, X[di], Y[di], Z[di], True)
    for i in range(1, di):
        u = i / di
        X[i], Y[i], Z[i] = frame1(u)

    """
    The final transition starts at k = npos // 2 + 2 * delay - di, with (X[k], Y[k], Z[k]),
    and ends at the next observing attitude (X1, Y1, Z1).
    """
    k = npos // 2 + 2 * delay - di
    q2a, frame2 = make_q01_and_interp_frame(X[k], Y[k], Z[k], X1, Y1, Z1, True)
    for i in range(1, di + 1):
        u = i / di
        X[k + i], Y[k + i], Z[k + i] = frame2(u)

    return X[1:], Y[1:], Z[1:]


# -----------------------------------------------------------------------------------------
def wrap_pi(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def get_obs_attitude(toSun, obj, roll=True):
    """
    Compute observing attitudes from the sequence of Sun directions toSun.
    Apply roll correction when requested.
    Before correction, Y is perpendicular to the Sun; pa is the uncorrected position angle.
    Parameters
    toSun : ndarray, shape=(N,3)
        Sun-direction unit vectors at each sample time.
    obj : dict
        Dictionary containing toObj, North, and target_pa.
    roll : bool

    Returns
    -------
    SatZ, SatY, SatX : ndarray, shape=(N, 3); ranmid : midpoint roll [deg]
    """
    n = len(toSun)
    SatZ = np.tile(obj["toObj"], (n, 1))     # Target direction
    SatY = normalize(np.cross(SatZ, toSun))  # Y axis perpendicular to the Sun direction
    SatX = normalize(np.cross(SatY, SatZ))

    pa = np.arctan2(SatX @ obj["North"], SatY @ obj["North"],)  # Position Angle

    if roll:
        target = np.radians(obj["target_pa"])
        desired_pa = np.where((-np.pi / 2 < pa) & (pa < np.pi / 2), target, np.pi + target,)

        ran = wrap_pi(desired_pa - pa)
        ranmid = np.degrees(ran[(n - 1) // 2])
        SatXr = np.empty_like(SatX)
        SatYr = np.empty_like(SatY)

        for i in range(n):
            q = quaternion.from_rotation_vector(SatZ[i] * ran[i])
            R = quaternion.as_rotation_matrix(q)

            SatXr[i] = normalize(R @ SatX[i])
            SatYr[i] = normalize(R @ SatY[i])

        SatX = SatXr
        SatY = SatYr
    else:
        ranmid = 0.

    return SatZ, SatY, SatX, ranmid


# -----------------------------------------------------------------------------------------
def center_jd(jd, satellite, T_orb, toObj):
    """
    Return the time of maximum target elevation within the orbit containing jd.
    """
    jdc1 = center_jd_step(jd, satellite, toObj, -T_orb / 2, T_orb / 2, 300)
    jdc2 = center_jd_step(jdc1, satellite, toObj, -300., 300., 30)
    jdc = center_jd_step(jdc2, satellite, toObj, -30., 30., 1)
    return jdc


# -----------------------------------------------------------------------------------------
def center_jd_step(jd, satellite, toObj, t1, t2, step):
    """
    Search from jd + t1 to jd + t2 in step-second increments for maximum target elevation.
    """
    n = int((t2 - t1) // step) + 1
    TargetElevations = np.empty(n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        warnings.simplefilter("ignore", AstropyWarning)
        for i in range(n):
            t = jd + timedelta(seconds=i * step + t1)
            TargetElevations[i] = TargetElevation(t, satellite, toObj)
        i = TargetElevations.argmax()
        jdc = Time(jd + timedelta(seconds=i * step + t1))

    # print("Sub",jdc.iso,ic,TargetElevations.max())
    return jdc


# -----------------------------------------------------------------------------------------
def TargetElevation(t, satellite, toObj):
    """
    Target elevation as seen from the satellite.
    """
    r, v = get_satellite_position(satellite, t)
    toSat = normalize(r)
    elevation = 90 - np.degrees(np.arccos(np.clip(np.dot(toObj, toSat), -1., 1.)))
    return elevation


# -----------------------------------------------------------------------------------------
def get_satellite_position(satellite, jd):
    """
    Compute satellite position and velocity.

    Parameters:
    - satellite: sgp4.model.Satrec
        SGP4 satellite record.
    - jd: astropy.Time or float
        Julian date.

    Returns:
    - position: ndarray
    - velocity: ndarray
    """
    if isinstance(jd, float):
        t = jd
    elif isinstance(jd, Time):
        t = jd.jd
    else:
        print('Type error')
        sys.exit()

    jdint = int(t)     # Integer part
    fr = t - jdint     # Fractional part
    # Compute position and velocity with SGP4.
    # Output coordinates are in the SGP4 TEME frame.
    e, position, velocity = satellite.sgp4(jdint, fr)
    if e != 0:
        raise RuntimeError(f"SGP4 error code: {e}")

    return position, velocity


# -----------------------------------------------------------------------------------------
def append_toobj(objs):
    """
    Add target direction and a nearby northward reference vector to each obj.
    """
    for obj in objs:
        target = obj["target"]
        ra, dec = target["ra"] * units.deg, target["dec"] * units.deg
        decN = (target["dec"] + 1) * units.deg           # One degree north
        objcoord = SkyCoord(ra=ra, dec=dec, frame='icrs')
        objdir = objcoord.cartesian.xyz.value
        toObj = normalize(objdir)
        objcoordN = SkyCoord(ra=ra, dec=decN, frame='icrs')
        objdirN = objcoordN.cartesian.xyz.value
        north = normalize(objdirN)
        obj['toObj'] = toObj
        obj['North'] = north


# -----------------------------------------------------------------------------------------
def get_jd_from_sun_lon(year, target_lon):
    """
    Find the Julian date corresponding to a solar ecliptic longitude in a given year.
    Parameters:
    - year: int
    - target_lon: astropy.u.degree

    Returns:
    - jd: astropy.Time
    """
    tlon = target_lon.value
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        warnings.simplefilter("ignore", AstropyWarning)
        start_date = Time(f"{year}-03-21") + timedelta(days=tlon * 365 / 360 - 20)
        end_date = Time(f"{year}-03-21") + timedelta(days=tlon * 365 / 360 + 20)
        result = minimize_scalar(lambda x: diff_sun_lon(x, target_lon),
                                 bounds=(start_date.jd, end_date.jd),
                                 options={'maxiter': 25}, method='bounded')
        jd = Time(result.x, format='jd')
    return jd


# -----------------------------------------------------------------------------------------
def diff_sun_lon(jd, target_lon):
    """
    Absolute difference between solar ecliptic longitude at jd and target_lon.

    Parameters:
    - jd: astropy.Time
    - jd,target_lon:  astropy.units.degree
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=ErfaWarning)
        warnings.filterwarnings('ignore', category=AstropyWarning)
        time = Time(jd, format='jd')
        sun = get_sun(time).transform_to(GeocentricTrueEcliptic())
    lon_diff = (sun.lon - target_lon).wrap_at(360 * units.deg)
    return abs(lon_diff.value)


# -----------------------------------------------------------------------------------------
def append_satrec(orbs):
    """
    Add an SGP4 Satrec and orbital period to each orbit dictionary.
    """
    for o in orbs:
        tle_line1, tle_line2, inclination_deg, mean_motion, T_orb = set_TLE(o)
        satellite = Satrec.twoline2rv(tle_line1, tle_line2)
        o["satellite"] = satellite
        o["Torb"] = T_orb


# -----------------------------------------------------------------------------------------
def set_TLE(orbit):
    """
    Construct a two-line element set (TLE).

    Parameters:
    - orbit: dict
        Dictionary containing the following fields:
        - epoch: str
            Reference epoch in UTC, formatted as "YYYY-MM-DD hh:mm:ss".
        - altitude: float
            Orbital altitude [km].
        - ltan: float
            Local time of ascending node [hours].
        - InternationalDesignator: str
            International designator: launch year, launch number, and piece.

    Returns:
    - tuple
        Tuple containing:

         - tle_line1: str
             First TLE line.
         - tle_line2: str
             Second TLE line.
         - inclination_deg : float
             Orbital inclination [deg].
         - mean_motion: float
              Mean motion [revolutions/day].
         - T_orb: float
               Orbital period [s].
    """
    # Constants
    mu = wgs84.mu  # Earth gravitational parameter [km^3/s^2]
    Re = wgs84.radiusearthkm  # Earth radius [km]
    J2 = wgs84.j2  # Earth J2 coefficient
    # J2 accounts for the oblateness correction to the gravitational potential
    # Days per tropical year
    days_per_year = 365.24219040
    seconds_per_day = 24 * 60 * 60
    # Orbital radius
    r = Re + orbit['altitude']
    # Compute inclination accounting for J2 nodal precession
    i_cos = (-(2 * np.pi / (days_per_year * seconds_per_day))
             / (1.5 * J2 * (Re / r)**2 * np.sqrt(mu / r**3)))
    if not (-1 <= i_cos <= 1):
        print('inclination angle error, i_cos ', i_cos)
        exit()
    inclination_rad = np.arccos(i_cos)  # Orbital inclination [rad]
    inclination_deg = np.degrees(inclination_rad)
    # Compute orbital period using Kepler's third law
    T_orb = 2 * np.pi * np.sqrt(r**3 / mu)  # [s]
    # Compute mean motion [revolutions/day]
    mean_motion = seconds_per_day / T_orb  # [revolutions/day]
    # Eccentricity (nearly circular orbit)
    eccentricity = 0.0000001  # Small nonzero eccentricity to avoid potential zero-eccentricity issues
    # Format eccentricity as seven unsigned digits
    eccentricity_str = f"{eccentricity:.7f}".split(".")[1].zfill(7)
    # Right ascension of the ascending node (RAAN)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        epoch = Time(orbit['epoch'], format='iso', scale='utc')
        raan = ltan_to_raan(orbit['ltan'], epoch)
        # Reference epoch
        tle_epoch = time_to_tle_epoch(epoch)

    # Construct TLE lines (68 characters without checksum)
    tle_line1 = (
        f"1 00000U {orbit['InternationalDesignator']:<8} "
        f"{tle_epoch}  .00000000  00000-0  00000-0 0  999"
    )
    tle_line2 = (
        f"2 00000 {inclination_deg:8.4f} {raan:8.4f} "
        f"{eccentricity_str} 000.0000 090.0000 "
        f"{mean_motion:11.8f}00000"
    )

    # Validate line lengths and append checksums
    if len(tle_line1) != 68 or len(tle_line2) != 68:
        raise ValueError("Invalid TLE line length")

    tle_line1 += tle_checksum(tle_line1)
    tle_line2 += tle_checksum(tle_line2)

    return tle_line1, tle_line2, inclination_deg, mean_motion, T_orb


# -----------------------------------------------------------------------------------------
def tle_checksum(line):
    """
    Calculate the TLE checksum for the first 68 characters.

    Digits contribute their numeric values, minus signs contribute 1,
    and all other characters contribute 0.
    """
    checksum = sum(int(c) for c in line if c.isdigit())
    checksum += line.count('-')
    return str(checksum % 10)


# -----------------------------------------------------------------------------------------
def ltan_to_raan(ltan, epoch):
    """
    Compute right ascension of ascending node (RAAN) from LTAN and epoch.

    Parameters:
    - ltan: float
        Local time of ascending node [hours].
    - epoch: astropy.time.Time

    Returns:
    - raan : float
        Right ascension of ascending node [deg].
    """
    # Apparent right ascension of the Sun
    solar_ra = get_sun(epoch).ra.to(units.deg).value
    # Convert LTAN to RAAN
    raan = (solar_ra + (ltan - 12.0) * 15) % 360
    return raan


# -----------------------------------------------------------------------------------------
def time_to_tle_epoch(time_obj):
    """
    Convert an astropy.time.Time object to TLE epoch format (YYDDD.DDDDDDDD).

    Parameters:
    - time_obj: astropy.time.Time
        The input time in UTC.

    Returns:
    - tle_epoch: str
        The TLE epoch string in YYDDD.DDDDDDDD format.
    """
    # Convert to Julian date to calculate day of the year
    year = time_obj.datetime.year
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        warnings.simplefilter("ignore", AstropyWarning)
        start_of_year = Time(f"{year}-01-01T00:00:00.000", format='isot', scale='utc')
        day_of_year = (time_obj - start_of_year).jd + 1  # Julian day difference + 1 for DOY

    # Extract the year in two-digit format
    year_short = int(str(year)[-2:])

    # Combine into TLE epoch format
    tle_epoch = f"{year_short:02d}{day_of_year:012.8f}"[0:14]
    return tle_epoch
