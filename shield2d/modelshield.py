#!/usr/bin/env python3
"""Two-dimensional ray tracing for a two-stage thermal shield.

Usage:
    modelshield.py [options]

Options:
    -a <rap>             Aperture radius [default: 300]
    -b <rbt>             Bottom radius [default: 199]
    -r <r1>              Intermediate radius [default: 261]
    -l <h1>              Intermediate height [default: 196]
    --height <h>         Total shield height [default: 600]
    --angle <degrees>    Reference elevation angle [default: 45]
    -o <figfile>         Save figure (PDF, PNG, etc.)
    --show               Display figure interactively
    -h --help            Show this help message
"""
import matplotlib.pyplot as plt
import numpy as np
from lib_shield import trace_ray
from docopt import docopt


def main(args):
    outfile = args['-o']
    rap = float(args['-a'])
    rbt = float(args['-b'])
    r1  = float(args['-r'])
    h1  = float(args['-l'])

    # Two-stage shield geometry is specified by the aperture, intermediate,
    # and bottom radii, together with the intermediate and total heights.
    # The reference elevation angle is used to highlight a traced ray.
    h = float(args["--height"])
    Tx = float(args["--angle"])
    if not (0 < rbt < r1 < rap and 0 < h1 < h):
        raise ValueError("Require 0 < rbt < r1 < rap and 0 < h1 < height")
    if not (0 < Tx < 90):
        raise ValueError("--angle must be between 0 and 90 degrees")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    plt.subplots_adjust(bottom=0.2)
    ax.set_aspect('equal')

    model = {"nlines": 6, "start_id": 0, "aperture_id": 5, "bottom_id": 2,
         0: {"A": np.array([-rap, h]), "B": np.array([-r1, h1])},
         1: {"A": np.array([-r1, h1]), "B": np.array([-rbt, 0.])},
         2: {"A": np.array([-rbt, 0.]), "B": np.array([rbt, 0])},
         3: {"A": np.array([rbt, 0]), "B": np.array([r1, h1])},
         4: {"A": np.array([r1, h1]), "B": np.array([rap, h])},
         5: {"A": np.array([rap, h]), "B": np.array([-rap, h])}}

    x = [model[0]["A"][0]]
    y = [model[0]["A"][1]]
    for i in range(model["nlines"]):
        x.append(model[i]["B"][0])
        y.append(model[i]["B"][1])

    Tok = 0
    for T in np.arange(20, 70.1, 1):
        P = model[0]["A"]
        U = np.array([np.cos(np.radians(T)), -np.sin(np.radians(T))])
        Ps, idy = trace_ray(model, P, U)

        if np.abs(T-Tx) < 0.1:
            ax.plot(Ps[:, 0], Ps[:, 1], color='blue', linestyle='solid')
            if idy == model["aperture_id"]:
                Tok = T
        elif idy == model["aperture_id"]:
            ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='solid')
            Tok = T
            print(f"T= {T} OK")
        else:
            ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='dotted')
    print(f'Max Tok = {Tok}')

    ax.plot(x, y, color='green')


    ax.set_title(f"Shield double type for elevation angle {Tx} degree")
    ax.text(0, -0.13,
            f"rap:{rap:.0f} rbt:{rbt:.0f} r1:{r1:.0f} h1:{h1:.0f} h:{h:.0f} ratio {rap/rbt:.2f}",
            transform=ax.transAxes)

    if outfile:
        from pathlib import Path
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile)

    if args["--show"] or not outfile:
        plt.show()
    plt.close(fig)


# -----------------------------------------------------------------------------------------
def rays_single(T, h, ax):
    Tr = np.radians(T)   # Convert the elevation angle to radians.
    # Geometric constraint: tan(T) = h / (r_ap + r_bt).
    rsum = h / np.tan(Tr)   # Sum of aperture and bottom radii.
    print(f"rsum {rsum}")
    # Choose the wall slope so a ray striking its lower edge reflects horizontally.
    # tan(T/2) = (r_ap - r_bt)/h
    rdif = np.tan(Tr / 2) * h
    print(f"rdif {rdif}")
    r_ap = (rdif + rsum) / 2
    print(f"r_ap {r_ap}")
    r_bt = rsum - r_ap
    print(f"r_bt {r_bt}")

    if ax is not None:
        model = {"nlines": 4, "start_id": 0, "aperture_id": 3, "bottom_id": 1,
                 0: {"A": np.array([-r_ap, h]), "B": np.array([-r_bt, 0.])},
                 1: {"A": np.array([-r_bt, 0.]), "B": np.array([r_bt, 0.])},
                 2: {"A": np.array([r_bt, 0.]), "B": np.array([r_ap, h])},
                 3: {"A": np.array([r_ap, h]), "B": np.array([-r_ap, h])}}

        x = [model[0]["A"][0]]
        y = [model[0]["A"][1]]
        for i in range(model["nlines"]):
            x.append(model[i]["B"][0])
            y.append(model[i]["B"][1])

        ax.plot(x, y, color='black')

        for Tt in np.arange(T - 1.01, T + 1.01, 1):   # Trace rays near the limiting elevation angle.
            P = model[0]["A"]
            U = np.array([np.cos(np.radians(Tt)), -np.sin(np.radians(Tt))])
            Ps, idy = trace_ray(model, P, U)
            if idy == model["aperture_id"]:
                ax.plot(Ps[:, 0], Ps[:, 1], color='blue', linestyle='solid')
            else:
                ax.plot(Ps[:, 0], Ps[:, 1], color='blue', linestyle='dotted')


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)
