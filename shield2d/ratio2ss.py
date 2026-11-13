#!/usr/bin/env python3
"""Plot aperture-to-bottom radius ratios for two-stage thermal shields.

Usage:
    ratio2ss.py [--height <h>] [--angles <angles>] [-o <file>] [--show]
    ratio2ss.py (-h | --help)

Options:
    --height <h>       Shield height [default: 600].
    --angles <angles>  Comma-separated limiting elevation angles in degrees
                       [default: 60,65,70].
    -o <file>          Output figure (PDF, PNG, etc.) [default: ratio2ss.pdf].
    --show             Display the figure interactively after saving.
    -h --help          Show this help message.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from docopt import docopt

from lib_shield import bestrbt


def parse_options():
    args = docopt(__doc__)
    try:
        height = float(args['--height'])
        angles = [float(x.strip()) for x in args['--angles'].split(',')]
        if height <= 0:
            raise ValueError('height must be positive')
        if not angles or any(not np.isfinite(x) or not (0 < x < 90) for x in angles):
            raise ValueError('angles must be between 0 and 90 degrees')
    except ValueError as exc:
        sys.exit(f'Invalid argument: {exc}')
    return height, angles, args['-o'], args['--show']


def main():
    height, angles, outfile, show = parse_options()
    fig, ax = plt.subplots()

    for angle in angles:
        rapbest, rbtbest, r1best, h1best, tdbest, raps, rbts = bestrbt(angle, height)
        ratio = rapbest / rbtbest
        ax.plot(rbts, np.asarray(raps) / rbts, label=rf'$\theta$={angle:g}')
        ax.plot([rbtbest], [ratio], marker='o', color='black')
        ax.text(rbtbest - 20, ratio - 0.15,
                f'rap:{rapbest:.1f} rbt:{rbtbest:.0f} ratio:{ratio:.1f}')

    ax.set_ylim(2, 5)
    ax.set_xlabel(r'$r_{bt}$')
    ax.set_ylabel(r'$r_{ap}/r_{bt}$')
    ax.set_title(f'Aperture-to-Bottom Ratio of a Two-Stage Shield (h={height:g})')
    ax.legend()
    fig.savefig(outfile, bbox_inches='tight')
    print(f'Saved: {outfile}')
    if show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
