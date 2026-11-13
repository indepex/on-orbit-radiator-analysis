#!/usr/bin/env python3
"""Plot Earth coupling over an orbit using attitude and outward-emission data.

Usage:
    earth_sat.py <attitude.npz> <outward.npz> [-o <plotfile>] [--show]
    earth_sat.py (-h | --help)

Options:
    -o <plotfile>  Save plot to this file (e.g. PNG or PDF).
    --show         Display the plot interactively, even when -o is given.
    -h --help      Show this help message.

The plotted coupling is the directional emission weight intercepted by the
Earth's apparent disk, divided by total emitted weight. It is not heat in W.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from docopt import docopt

from lib_sat import earth_couple, load_outward


def main(args):
    attitude_file = args['<attitude.npz>']
    outward_file = args['<outward.npz>']
    plot_file = args['-o']

    (modelfile, ntot, total_emission, absorbed, dirs, dir_emission, npt, area
     ) = load_outward(outward_file)
    dirs = np.asarray(dirs)
    dir_emission = np.asarray(dir_emission)
    if dirs.ndim != 2 or dirs.shape[1] != 3 or dir_emission.shape != (len(dirs),):
        raise ValueError('Inconsistent Dirs and dir_emission shapes in outward NPZ')
    if not np.isfinite(total_emission) or total_emission <= 0:
        raise ValueError('total_emission must be positive and finite')

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7']
    linestyles = ['-', '--', ':', '-.']
    ax.set_prop_cycle(cycler(linestyle=linestyles) * cycler(color=colors))

    maxima = []
    altitudes = []
    ltans = []
    with np.load(attitude_file) as data:
        ncase = int(data['ncase'])
        if ncase < 1:
            raise ValueError('Attitude NPZ contains no cases')

        for k in range(ncase):
            p = f'case{k:03d}'
            altitude = float(data[f'{p}_altitude'])
            ltan = float(data[f'{p}_ltan'])
            roll = float(data[f'{p}_roll_angle'])
            slon = float(data[f'{p}_sun_longitude'])
            theta = np.asarray(data[f'{p}_theta'])
            to_sat = np.asarray(data[f'{p}_toSat'])
            sat_x = np.asarray(data[f'{p}_SatX'])
            sat_y = np.asarray(data[f'{p}_SatY'])
            sat_z = np.asarray(data[f'{p}_SatZ'])
            n = len(theta)
            if theta.ndim != 1 or any(v.shape != (n, 3) for v in (to_sat, sat_x, sat_y, sat_z)):
                raise ValueError(f'Inconsistent attitude arrays for {p}')

            earth_cpl, space = earth_couple(
                dirs, dir_emission, to_sat, sat_x, sat_y, sat_z, altitude
            )
            coupling = earth_cpl / total_emission
            mean_cpl = float(np.mean(coupling))
            maxima.append(float(np.max(coupling)))
            altitudes.append(altitude)
            ltans.append(ltan)
            ax.plot(theta, coupling,
                    label=f'Mean: {mean_cpl:.3f}  Slon: {slon:.1f}  Roll: {roll:.1f}')
            print(f'{k:03d} Mean {mean_cpl:.3f} Slon:{slon:5.1f} '
                  f'LTAN:{ltan:4.3f} Altitude:{altitude:.0f} Roll:{roll:.1f}')

    title = f'Earth coupling variation — {Path(str(modelfile)).name}'
    if np.allclose(altitudes, altitudes[0]) and np.allclose(ltans, ltans[0]):
        title += f'  Alt {altitudes[0]:.0f} km  LTAN {ltans[0]:.2f}'
    title += f'  Area {float(area) / 1e6:.3g} m²'
    ax.set_title(title)
    ax.set_xlabel('Satellite phase angle [degree]')
    ax.set_ylabel('Earth coupling fraction')
    ax.set_xlim(0, 360)
    ax.set_ylim(bottom=0, top=max(0.01, max(maxima) * 1.08))
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize='small')
    fig.tight_layout()

    if plot_file:
        output = Path(plot_file)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=180, bbox_inches='tight')
        print(f'Saved: {output}')
    if args['--show'] or not plot_file:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main(docopt(__doc__))
