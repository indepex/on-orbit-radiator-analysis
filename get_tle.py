#!/usr/bin/env python3
"""
get TLE

    Usage:
        get_tle.py --orb <orbit_file>

    Options:
        --orb <orbit_file>  orbit specifications
"""
from docopt import docopt
from lib_sat import (set_TLE)
from lib_utils import (load_json5)


def main(args):
    orb = load_json5(args['--orb'], "Orbit")
    tle_line1, tle_line2, inclination_deg, mean_motion, T_orb = set_TLE(orb)
    print(tle_line1)
    print(tle_line2)


if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)
