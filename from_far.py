#!/usr/bin/env python3
r"""
Trace rays incident on a three-dimensional spacecraft model.

Usage:
    from_far.py <level> <model.json5> [-t <theta>] [-p <phi>] [-c] [-o <plotfile>] [--plot-step <n>]
    from_far.py (-h | --help)

Arguments:
    <level>          0: display the model only
                     1: register and inspect surfaces in the C library
                     2: trace one incident ray
                     3: trace a 41 x 41 grid of incident rays
                     4: trace a 35 x 51 grid in a selected region
    <model.json5>    Spacecraft model parameters

Options:
    -t <theta>       Incident zenith angle from +Z, in degrees [default: 0]
    -p <phi>         Incident azimuth from +X toward +Y, in degrees [default: 0]
    --plot-step <n>  Plot every nth ray in both grid directions [default: 1]
    -c               Use the predefined camera position
    -o <plotfile>    Save the rendered plot
    -h --help        Show this help
"""

import numpy as np
import pyvista as pv
from docopt import docopt

from lib_model import load_payload_model, generate_block, id2name, flagname, name2id
from lib_radcpl import set_surface, gen_rays
from lib_view import view_block
from lib_utils import normalize


def incident_direction(theta, phi):
    """Return the propagation direction of a ray arriving from (theta, phi)."""
    th, ph = np.radians([theta, phi])
    return -np.array([np.sin(th) * np.cos(ph),
                      np.sin(th) * np.sin(ph),
                      np.cos(th)])


def trace_rays(model, level, theta, phi):
    """Trace incident rays, preserving the original launch geometries."""
    set_surface(model)
    sid = 0  # Dummy surface for rays launched outside the spacecraft
    z_axis = np.array([0, 0, 1])
    direction = incident_direction(theta, phi)
    distance = 3000

    if level == 2:
        starts = -distance * direction + 1303 * z_axis
        return gen_rays(starts, sid, direction, id_out=True, shell=True)

    center = -distance * direction + 1003 * z_axis
    local_x = normalize(np.cross(direction, z_axis))
    local_y = normalize(np.cross(direction, local_x))

    if level == 3:
        starts = [center + ix * 50 * local_x + iy * 50 * local_y
                  for ix in range(-20, 21) for iy in range(-20, 21)]
    else:  # level == 4: fixed grid for detailed design inspection
        starts = [center + (ix + 10) * 20 * local_x
                  + (iy - 8) * 20 * local_y
                  for ix in range(-5, 30) for iy in range(-50, 1)]

    return gen_rays(starts, sid, direction, id_out=True,
                    shell=True, shell_rad=4500)


def ray_blocks(model, radiator, points, ids, flags, level, plot_step):
    """Create ray line blocks and print the original hit summary."""
    blocks, colors, opacities = [], [], []
    direct, radiator_hits, other_hits = 0, 0, 0

    # Number of grid points along the inner-loop (Y) direction.
    grid_width = {3: 41, 4: 51}.get(level)

    for k, ray in enumerate(points):
        status = flagname(flags[k])
        absorbed = status == 'HIT_ABSORB'
        last_id = ids[k][-1] if absorbed else None

        # Keep the full-ray statistics independent of the display sampling.
        plot_ray = (grid_width is None
                    or ((k // grid_width) % plot_step == 0
                        and (k % grid_width) % plot_step == 0))

        if absorbed:
            print(f'absorbed on {id2name(model, last_id)}')

        if plot_ray:
            for i in range(len(ray) - 1):
                blocks.append(pv.Line(ray[i], ray[i + 1]))
                if absorbed:
                    colors.append('red' if last_id == radiator else 'orange')
                elif len(ray) == 2:
                    colors.append('green')
                else:
                    colors.append('blue')
                opacities.append(1)

        if absorbed:
            if last_id == radiator:
                radiator_hits += 1
            else:
                other_hits += 1
        if len(ray) == 2 and status == 'HIT_NONE':
            direct += 1

    print(f'Direct {direct} Radiator {radiator_hits} '
          f'Absorbed {other_hits} Total {len(points)}')
    return blocks, colors, opacities


def main(modelfile, level, savefile, cam, van, psc, theta, phi, plot_step):
    model, cubes = load_payload_model(modelfile)
    cblk, ccol, copa = generate_block(cubes)
    mblk, mcol, mopa = generate_block(model)
    rblk, rcol, ropa = [], [], []

    if level == 1:
        set_surface(model, debug=True)
        print('Push data into C-library')
    elif level >= 2:
        radiator = name2id(model, 'Radiator')
        points, ids, flags, final_dirs = trace_rays(model, level, theta, phi)
        rblk, rcol, ropa = ray_blocks(model, radiator, points, ids, flags, level, plot_step)

    combined = pv.MultiBlock(cblk + mblk + rblk)
    view_block(combined, ccol + mcol + rcol, copa + mopa + ropa,
               savefile, campos=cam, viewan=van, pscale=psc)


def parse_arguments():
    args = docopt(__doc__)
    level = int(args['<level>'])
    if level not in range(5):
        raise SystemExit('Error: <level> must be an integer from 0 to 4.')

    plot_step = int(args['--plot-step'])
    if plot_step < 1:
        raise SystemExit('Error: --plot-step must be a positive integer.')

    if args['-c']:
        cam = [(-9098.329043666137, 12218.910665226284, 11449.936043069698),
               (-341.9481946375063, -43.9962589282386, 869.0264465392593),
               (0.06690868765602145, -0.6239943938838902, 0.7785590689970966)]
        van, psc = 30., 2500.
    else:
        cam, van, psc = None, None, None

    return (args['<model.json5>'], level, args['-o'], cam, van, psc,
            float(args['-t']), float(args['-p']), plot_step)


if __name__ == '__main__':
    main(*parse_arguments())
