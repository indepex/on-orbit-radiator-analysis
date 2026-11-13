#!/usr/bin/env python3
r"""
Calculate the directional radiation coupling of the spacecraft radiator.

Usage:
    outward_emission.py <model.json5> <num> -o <outfile> [--grid-step <step>]
    outward_emission.py (-h | --help)

Arguments:
    <model.json5>       Spacecraft model definition.
    <num>               Number of emission directions per starting point.

Options:
    -o <outfile>        Output NPZ file.
    --grid-step <step>  Radiator sampling grid spacing [default: 10].
    -h --help           Show this help.
"""

from pathlib import Path

import numpy as np
from docopt import docopt

from lib_model import load_payload_model, name2id
from lib_radcpl import set_surface, fibonacci_sphere_half, count_out_dirs


def parse_args():
    args = docopt(__doc__)
    num = int(args['<num>'])
    step = float(args['--grid-step'])
    if num < 1:
        raise SystemExit('Error: <num> must be a positive integer.')
    if not np.isfinite(step) or step <= 0:
        raise SystemExit('Error: --grid-step must be a positive finite number.')
    return args['<model.json5>'], num, args['-o'], step


def main(modelfile, num, savefile, dr):
    model, cubes = load_payload_model(modelfile)
    set_surface(model)

    # Use the same hemisphere directions for launching and binning rays.
    directions = fibonacci_sphere_half(num, True)

    sid = name2id(model, 'Radiator')
    bid = name2id(model, 'TelBaffle')
    bx, by = model[bid]['p0'][0], model[bid]['p0'][1]
    br = model[bid]['radius']

    # Radiator is a rectangle in the XY plane; exclude the baffle footprint.
    x0, x1 = sorted((model[sid]['p0'][0], model[sid]['p1'][0]))
    y0, y1 = sorted((model[sid]['p0'][1], model[sid]['p2'][1]))
    xs = np.arange(x0 + dr / 2, x1, dr)
    ys = np.arange(y0 + dr / 2, y1, dr)

    starts = []
    area = 0.0
    for y in ys:
        for x in xs:
            if (x - bx) ** 2 + (y - by) ** 2 > br ** 2:
                starts.append(np.array([x, y, model[sid]['p0'][2]]))
                area += dr * dr

    npoint = len(starts)
    if npoint == 0:
        raise SystemExit('Error: no radiator sampling points were generated.')
    total = len(directions) * npoint
    print(f'Start points: {npoint}; directions: {len(directions)}; rays: {total}')

    dir_emission, total_emission, absorbed_emission = count_out_dirs(
        starts, sid, directions
    )
    print(f'Absorbed {absorbed_emission} Total {total_emission}')

    output = Path(savefile)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output,
        modelfile=modelfile,
        num_rays_total=total,
        total_emission=total_emission,
        absorbed_emission=absorbed_emission,
        Dirs=directions,
        dir_emission=dir_emission,
        Npoint=npoint,
        Area=area,
    )


if __name__ == '__main__':
    main(*parse_args())
