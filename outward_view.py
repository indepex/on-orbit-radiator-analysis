#!/usr/bin/env python3
r"""
Plot the directional radiation coupling calculated by outward_emission.py.

Usage:
    outward_view.py <datafile> [-o <plotfile>] [-e <el>] [-a <az>] [--show]
    outward_view.py (-h | --help)

Arguments:
    <datafile>       Input NPZ file from outward_emission.py.

Options:
    -o <plotfile>    Save the plot to a file (e.g. PNG or PDF).
    -e <el>          Viewing elevation angle in degrees.
    -a <az>          Viewing azimuth angle in degrees.
    --show           Display an interactive plot.
    -h --help       Show this help.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from docopt import docopt

from lib_sat import load_outward


def parse_args():
    args = docopt(__doc__)
    if not args['-o'] and not args['--show']:
        raise SystemExit('Error: specify -o <plotfile> and/or --show.')
    el = float(args['-e']) if args['-e'] is not None else None
    az = float(args['-a']) if args['-a'] is not None else None
    return args['<datafile>'], args['-o'], el, az, args['--show']


def main(datafile, plotfile, el, az, show):
    (modelfile, ntot, total_emission, absorbed, dirs, dir_emission, npt, area
     ) = load_outward(datafile)

    peak = np.max(dir_emission)
    if peak <= 0:
        raise SystemExit('Error: no positive directional emission values found.')
    relative_emission = dir_emission / peak

    fig = plt.figure(figsize=(8.0, 6.5))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    scatter = ax.scatter(
        dirs[:, 0], dirs[:, 1], dirs[:, 2],
        c=relative_emission, cmap='hot', vmin=0.0, vmax=1.0,
        marker='o', s=2, linewidths=0, rasterized=True,
    )
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.68, pad=0.12)
    colorbar.set_label('Relative directional emission (max = 1)')
    ratio = absorbed / total_emission * 100 if total_emission else 0.0
    ax.set_title('Directional Emission from the Radiator', pad=18)
    fig.suptitle(
        f'{Path(modelfile).name}  |  Total: {total_emission:.0f}  |  '
        f'Absorbed: {absorbed:.0f} ({ratio:.1f}%)',
        fontsize=10, y=0.96,
    )
    if el is not None and az is not None:
        ax.view_init(elev=el, azim=az)

    if plotfile:
        output = Path(plotfile)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=180, bbox_inches='tight')
        print(f'Saved {output}')
    if show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main(*parse_args())
