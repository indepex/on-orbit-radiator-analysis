#!/usr/bin/env python3
"""
Evaluate directional coupling to Earth radiation using a simplified dice satellite model.

    Usage:
        earth_coupling.py [options]

    Options:
        --orb <orbit_file>  Orbit configuration file [default: ../json/orbit.json5]
        --ltan <ltan>       Single LTAN value; defaults to the value in --orb
        --alt <alt>         Single altitude value; defaults to the value in --orb
        --obj <obj_file>    Observation target configuration file [default: ../json/GC_pos.json5]
        --sln <slon>        Single solar ecliptic longitude; defaults to the value in --obj
        --tpa <target_pa>   Target position-angle offset [default: -2] (spring: -2, autumn: 178)
        -o <plot_file>
        --dly <delay>      Observation-duration reduction parameter (default: 0; 1 is supported)


"""
from docopt import docopt
import numpy as np
from sgp4.earth_gravity import wgs84
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

from lib_utils import load_json5, normalize
from lib_sat import append_satrec, append_toobj, get_attitude
from lib_radcpl import fibonacci_sphere


def main(args):
    # Generate approximately uniformly distributed sampling directions.
    num = 200000
    direction = fibonacci_sphere(num)
    # Orbit and attitude calculation settings.
    npos = 72   # Number of orbital intervals (multiple of 4); produces npos + 1 samples.
    orb_d, obj_d = setup_data(args)            # Load orbit and target configuration dictionaries from CLI arguments.
    # The append_* functions accept lists of dictionaries; wrap each single case in a list.
    append_satrec([orb_d])
    append_toobj([obj_d])
    delay = int(args['--dly']) if args['--dly'] is not None else 0
    if args['--dly'] is not None:
        print(f"delay {delay}")

    alt, ltan, pa = orb_d['altitude'], orb_d['ltan'], obj_d['target_pa']
    slon = obj_d['sun_longitude']
    if 90 < slon < 270:
        pa = 180 + pa

    # get_attitude also returns time, position, and velocity.
    toSat, theta, SatX, SatY, SatZ, jdc, toSun, roll, *_ = get_attitude(
        orb_d, obj_d, npos, delay=delay
    )
    title = f"Altitude {alt:.0f}km ltan {ltan:.2f} Roll {roll:.2f} PA {pa:.0f} S.lon{slon}"
    print(title)
    counter = count_earth_direction(direction, SatX, SatY, SatZ, toSat, alt)
    plot_lonlat(direction, counter, title, args['-o'])


def plot_lonlat(direction, counter, title, outfile):
    # Longitude-latitude map with +Z as the north pole and -X at zero longitude.
    x, y, z = direction[:, 0], direction[:, 1], direction[:, 2]
    lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    lon = np.degrees(np.arctan2(y, -x))
    vmin, vmax = 0.0001, 11

    print(f"Counter max {np.max(counter)}")
    levels = np.linspace(vmin, vmax, 21)
    lon_grid = np.linspace(-180, 180, 361)
    lat_grid = np.linspace(-90, 90, 181)
    LON, LAT = np.meshgrid(lon_grid, lat_grid)
    C = griddata(points=(lon, lat), values=counter, xi=(LON, LAT), method="linear",)
    C = gaussian_filter(C, sigma=1.0)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)

    Cplot = np.where(C >= vmin, C, np.nan)
    ax.contour(LON, LAT, C, levels=levels,  vmin=vmin, vmax=vmax)
    # cs = ax.contour(LON, LAT, C, levels=levels,  vmin=vmin, vmax=vmax)
    # ax.clabel(cs, inline=True, fontsize=8)

    ax.contourf(LON, LAT, Cplot, levels=levels, vmin=vmin, vmax=vmax)
    # cs = ax.contourf(LON, LAT, Cplot, levels=levels, vmin=vmin, vmax=vmax)
    # cbar = fig.colorbar(cs, ax=ax, orientation="vertical")

    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 90)
    ax.set_xlabel("Longitude [deg], 0 deg = -SatX")
    ax.set_ylabel("Latitude [deg], north pole = +SatZ")

    ax.set_title(title)
    ax.axhline(y=0, color='r', linestyle=':')
    ax.grid(True, linestyle=":")

    if outfile is not None:
        fig.savefig(outfile, bbox_inches="tight")

    plt.show()


def count_earth_direction(direction, SatX, SatY, SatZ, toSat, alt):
    # Accumulate the positive Z-direction cosine for directions intersecting Earth.
    counter = np.zeros(len(direction), dtype=float)

    earth_radius_km = wgs84.radiusearthkm
    alpha = np.arcsin(earth_radius_km / (earth_radius_km + alt))
    cos_alpha = np.cos(alpha)

    # Test Earth visibility in spacecraft coordinates at each orbital sample.
    for i in range(len(SatX)):
        e_body = to_earth_in_sat_coord(SatX[i], SatY[i], SatZ[i], toSat[i])
        cosang = direction @ e_body
        mask = (cosang >= cos_alpha) & (direction[:, 2] > 0)
        counter[mask] += direction[mask, 2]

    return counter


def to_earth_in_sat_coord(SatX, SatY, SatZ, toSat):
    """
    Convert the inertial Earth direction (-toSat) to spacecraft-frame components.
    """
    earth_inertial = normalize(-toSat)

    return np.array([
        np.dot(SatX, earth_inertial),
        np.dot(SatY, earth_inertial),
        np.dot(SatZ, earth_inertial),
    ])


def single_float(arg, default, name):
    """Convert a CLI value or JSON5 default to a single floating-point value."""
    value = arg if arg is not None else default
    if isinstance(value, (list, tuple, np.ndarray)):
        if np.asarray(value).size != 1:
            raise ValueError(f"{name} must be a single value: {value}")
        value = np.asarray(value).item()
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a single numeric value: {value}") from exc


def setup_data(args):
    print(f"orbit file {args['--orb']}")
    orb = load_json5(args['--orb'], "Orbit")
    obj = load_json5(args['--obj'], "Object")

    orb['ltan'] = single_float(args['--ltan'], orb['ltan'], '--ltan')
    orb['altitude'] = single_float(args['--alt'], orb['altitude'], '--alt')
    obj['sun_longitude'] = single_float(args['--sln'], obj['sun_longitude'], '--sln')
    obj['target_pa'] = single_float(args['--tpa'], obj['target_pa'], '--tpa')
    return orb, obj


def add_attitude_result(data, k, orb, obj, toSat, theta, SatX, SatY, SatZ, jdc, toSun, roll_angle):
    p = f"case{k:03d}"

    # Orbit parameters
    data[f"{p}_altitude"] = np.array(orb["altitude"])
    data[f"{p}_ltan"] = np.array(orb["ltan"])
    data[f"{p}_Torb"] = np.array(orb["Torb"])
    # Observation target
    data[f"{p}_target_name"] = np.array(obj["target"].get("name", ""))
    data[f"{p}_toObj"] = np.asarray(obj["toObj"])
    # Observation conditions
    data[f"{p}_year"] = np.array(obj["year"])
    data[f"{p}_sun_longitude"] = np.array(obj["sun_longitude"])
    data[f"{p}_target_pa"] = np.array(obj["target_pa"])
    data[f"{p}_roll_angle"] = np.array(roll_angle)
    data[f"{p}_jdc_jd"] = np.array(jdc.jd)
    # Orbit-sampled data
    data[f"{p}_theta"] = np.asarray(theta)
    data[f"{p}_toSat"] = np.asarray(toSat)
    data[f"{p}_toSun"] = np.asarray(toSun)
    data[f"{p}_SatX"] = np.asarray(SatX)
    data[f"{p}_SatY"] = np.asarray(SatY)
    data[f"{p}_SatZ"] = np.asarray(SatZ)


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)
