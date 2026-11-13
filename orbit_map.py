#!/usr/bin/env python3
"""
Plot a satellite ground track from a TLE using SGP4.

Usage:
    orbit_map.py --tle <file> [--start <datetime>] [--days <days>] [--step <seconds>] [-o <outfile>]

Options:
    -h --help           Show this help message.
    --tle <file>        Input TLE file (two lines).
    --start <datetime>  Start time in UTC: "YYYY-MM-DD HH:MM:SS".
                        Default: TLE epoch.
    --days <days>       Calculation duration in days [default: 1].
    --step <seconds>    Time step in seconds [default: 60].
    -o <outfile>        Output figure (PNG, PDF, etc.).
                        If omitted, display interactively.
"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from docopt import docopt
from sgp4.api import Satrec
from astropy.time import Time, TimeDelta
from astropy.coordinates import TEME, ITRS, CartesianRepresentation
from astropy import units as u
from astropy.utils.exceptions import AstropyWarning
from erfa import ErfaWarning


# -----------------------------------------------------------------------------------------
def main(args):
    tle1, tle2 = read_tle(args["--tle"])
    satellite = Satrec.twoline2rv(tle1, tle2)

    if satellite.error != 0:
        raise RuntimeError(f"SGP4 initialization error: {satellite.error}")

    if args["--start"] is None:
        start = Time(
            satellite.jdsatepoch + satellite.jdsatepochF,
            format="jd",
            scale="utc",
        )
    else:
        start = Time(args["--start"], format="iso", scale="utc")

    days = float(args["--days"])
    step = float(args["--step"])

    if days <= 0 or step <= 0:
        raise ValueError("--days and --step must be positive")

    times, lat, lon = calculate_orbit(satellite, start, days, step)

    plot_orbit(times, lat, lon, satellite, days, args["-o"])


# -----------------------------------------------------------------------------------------
def read_tle(filename):
    """Read and validate two TLE lines."""

    with open(filename, "r", encoding="ascii") as f:
        lines = [line.rstrip("\r\n") for line in f if line.strip()]

    if len(lines) != 2:
        raise ValueError("TLE file must contain exactly two nonempty lines")

    tle1, tle2 = lines

    if not tle1.startswith("1 ") or not tle2.startswith("2 "):
        raise ValueError("Invalid TLE line numbers")

    for line in (tle1, tle2):
        if len(line) != 69:
            raise ValueError("TLE lines must be 69 characters long")

        checksum = sum(int(c) for c in line[:68] if c.isdigit())
        checksum += line[:68].count("-")

        if checksum % 10 != int(line[68]):
            raise ValueError("Invalid TLE checksum")

    return tle1, tle2


# -----------------------------------------------------------------------------------------
def calculate_orbit(satellite, start, days, step):
    """
    Propagate the orbit using SGP4 and convert positions to latitude/longitude.

    Returns
    -------
    times : astropy.time.Time
        Sample times in UTC.
    lat : ndarray
        Geodetic latitude in degrees.
    lon : ndarray
        Longitude in degrees.
    """

    duration = days * 86400.0
    offsets = np.arange(0.0, duration, step)

    if len(offsets) == 0 or offsets[-1] < duration:
        offsets = np.append(offsets, duration)

    times = start + TimeDelta(offsets, format="sec")

    # SGP4: TLE -> TEME position in km
    errors, positions, velocities = satellite.sgp4_array(
        times.jd1, times.jd2
    )

    if np.any(errors != 0):
        raise RuntimeError(
            f"SGP4 propagation errors: {np.unique(errors[errors != 0])}"
        )

    # TEME -> ITRS (Earth-fixed coordinates)
    xyz = CartesianRepresentation(
        positions[:, 0] * u.km,
        positions[:, 1] * u.km,
        positions[:, 2] * u.km,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        warnings.simplefilter("ignore", AstropyWarning)

        teme = TEME(xyz, obstime=times)
        itrs = teme.transform_to(ITRS(obstime=times))

    location = itrs.earth_location
    lat = np.asarray(location.lat.deg)
    lon = np.asarray(location.lon.deg)

    return times, lat, lon


# -----------------------------------------------------------------------------------------
def plot_orbit(times, lat, lon, satellite, days, outfile):
    """Plot a ground track with the first revolution highlighted."""

    fig = plt.figure(figsize=(12, 6.3))
    ax = fig.add_subplot(
        111,
        projection=ccrs.Mollweide(central_longitude=0),
    )

    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="0.95")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.gridlines(linewidth=0.4, alpha=0.5)

    # Orbital period derived from the TLE mean motion.
    mean_motion = satellite.no_kozai * 1440.0 / (2.0 * np.pi)
    period = 86400.0 / mean_motion

    elapsed = (times - times[0]).to_value("sec")
    first_orbit = elapsed <= period

    # Draw the complete ground track.
    plot_track(
        ax, lon, lat,
        color="tab:blue",
        linewidth=0.8,
        zorder=2,
    )

    # Highlight the first orbital revolution.
    plot_track(
        ax, lon[first_orbit], lat[first_orbit],
        color="tab:red",
        linewidth=1.8,
        zorder=3,
    )

    inclination = np.degrees(satellite.inclo)

    ax.set_title(
        "Ground Track of a Sun-Synchronous Orbit",
        fontsize=14,
        pad=14,
    )

    info = (
        f"Start: {times[0].utc.iso} UTC    "
        f"Duration: {days:.1f} days\n"
        f"Inclination: {inclination:.4f} deg    "
        f"Mean motion: {mean_motion:.8f} rev/day"
    )

    fig.text(
        0.5, 0.065, info,
        ha="center",
        va="center",
        fontsize=9,
    )

    fig.subplots_adjust(
        left=0.025,
        right=0.975,
        top=0.91,
        bottom=0.14,
    )

    if outfile:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
        print(f"Saved: {outfile}")
    else:
        plt.show()

    plt.close(fig)


# -----------------------------------------------------------------------------------------
def plot_track(ax, lon, lat, **kwargs):
    """
    Plot a ground track without connecting points across the
    +/-180 degree longitude boundary.
    """

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    # Split into continuous longitude segments.
    breaks = np.where(np.abs(np.diff(lon)) > 180.0)[0] + 1
    segments = np.split(np.arange(len(lon)), breaks)

    for segment in segments:
        if len(segment) < 2:
            continue

        ax.plot(
            lon[segment],
            lat[segment],
            transform=ccrs.PlateCarree(),
            **kwargs,
        )


# -----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main(docopt(__doc__))
