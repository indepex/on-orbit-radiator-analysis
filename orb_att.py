#!/usr/bin/env python3
"""
orbit and attitude

    Usage:
        orb_att.py -s <save_file> [options]

    Options:
        -s <save_file>      Output NPZ file for the calculation results
        --orb <orbit_file>  Orbit configuration file [default: ../json/orbit.json5]
        --ltan <ltan_list>  LTAN value or comma-separated list (default: from --orb)
        --alt <alt_list>    Altitude value or comma-separated list (default: from --orb)
        --obj <obj_file>    Target object configuration file [default: ../json/GC_pos.json5]
        --sln <slon_list>   Solar ecliptic longitude or comma-separated list (default: from --obj)
        --tpa <target_pa>   Target position-angle offset [default: -2] (spring: -2; autumn: 178)
        --dly <delay>      Parameter to shorten the observation interval (default: 0; 1 allowed)

"""
from docopt import docopt
import numpy as np
import sys
import shlex
from lib_utils import load_json5, expand_dict, list_rewrite
from lib_sat import (append_satrec, append_toobj, get_attitude)


def main(args):
    npos = 72   # Number of sampling intervals per orbit; must be divisible by 4. Arrays contain npos + 1 samples.
    orb_d, obj_d = setup_data(args)            # Load orbit and target configuration dictionaries from the arguments
    orbs = expand_dict(orb_d, ['ltan', 'altitude'])   # Expand into a list of dictionaries for all orbit parameter combinations
    append_satrec(orbs)  # Add an SGP4 Satrec object and orbital period to each orbit dictionary
    objs = expand_dict(obj_d, ['sun_longitude', 'target_pa'])  # Expand target and observation parameter combinations
    append_toobj(objs)  # Add toObj and North vectors to each target dictionary
    if args['--dly']:
        delay = int(args['--dly'])
        print(f"delay {delay}")
    else:
        delay = 0
    data = {}
    case_id = 0
    for orb in orbs:
        for obj in objs:
            (toSat, theta, SatX, SatY, SatZ, jdc, toSun, roll_angle,
             time_jd, position, velocity) = get_attitude(orb, obj, npos, delay=delay)

            print(orb)

            add_attitude_result(
                data, case_id, orb, obj, toSat, theta, SatX, SatY, SatZ, jdc, toSun, roll_angle,
                time_jd, position, velocity
            )
            case_id += 1

    data["ncase"] = np.array(case_id)
    data["npos"] = np.array(npos)
    # Save the results
    data["coordinate_system"] = np.array(
        "Position and velocity: TEME (SGP4). "
        "Sun direction: geocentric Astropy solar direction. "
        "Target direction: ICRS. "
        "Attitude axes: approximate common inertial frame; "
        "TEME, GCRS and ICRS differences are neglected."
    )
    data["command_line"] = np.array(shlex.join(sys.argv))
    np.savez_compressed(args["-s"], **data)


def setup_data(args):
    # Load configuration data, using file defaults unless overridden by CLI arguments
    # Orbit configuration
    print(f"orbit file {args['--orb']}")
    orb = load_json5(args['--orb'], "Orbit")
    # Allow LTAN to be overridden from the command line
    orb['ltan'] = list_rewrite(args['--ltan'], orb['ltan'])
    # Allow altitude to be overridden from the command line
    orb['altitude'] = list_rewrite(args['--alt'], orb['altitude'])
    # Target object configuration
    obj = load_json5(args['--obj'], "Object")
    # Allow solar ecliptic longitude to be overridden from the command line
    obj['sun_longitude'] = list_rewrite(args['--sln'], obj['sun_longitude'])
    # Allow target position angle to be overridden from the command line
    obj['target_pa'] = list_rewrite(args['--tpa'], obj['target_pa'])

    return orb, obj


def add_attitude_result(data, k, orb, obj, toSat, theta,
                        SatX, SatY, SatZ, jdc, toSun, roll_angle,
                        time_jd, position, velocity):
    """
    Append one orbit of orbit and attitude results to the NPZ output dictionary.

    time_jd : JD (UTC)
    position: TEME [km]
    velocity: TEME [km/s]

    The attitude axes are calculated by approximately treating TEME, GCRS,
    and ICRS as a common inertial reference frame.
    """
    p = f"case{k:03d}"

    # Orbit parameters
    data[f"{p}_altitude"] = np.array(orb["altitude"])
    data[f"{p}_ltan"] = np.array(orb["ltan"])
    data[f"{p}_Torb"] = np.array(orb["Torb"])

    # Target
    data[f"{p}_target_name"] = np.array(obj["target"].get("name", ""))
    data[f"{p}_toObj"] = np.asarray(obj["toObj"])

    # Observation parameters
    data[f"{p}_year"] = np.array(obj["year"])
    data[f"{p}_sun_longitude"] = np.array(obj["sun_longitude"])
    data[f"{p}_target_pa"] = np.array(obj["target_pa"])
    data[f"{p}_roll_angle"] = np.array(roll_angle)
    data[f"{p}_jdc_jd"] = np.array(jdc.jd)

    # Per-orbit sampled data
    data[f"{p}_theta"] = np.asarray(theta)
    data[f"{p}_toSat"] = np.asarray(toSat)
    data[f"{p}_toSun"] = np.asarray(toSun)
    data[f"{p}_SatX"] = np.asarray(SatX)
    data[f"{p}_SatY"] = np.asarray(SatY)
    data[f"{p}_SatZ"] = np.asarray(SatZ)
    data[f"{p}_time_jd"] = np.asarray(time_jd)
    data[f"{p}_position"] = np.asarray(position)
    data[f"{p}_velocity"] = np.asarray(velocity)


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)
