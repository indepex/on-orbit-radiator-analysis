#!/usr/bin/env python3
"""
Compare single- and two-stage thermal shield cross sections.

Usage:
    shield.py --angle <degrees> [-o <figfile>] [--height <height>]
    shield.py (-h | --help)

Options:
    --angle <degrees>  Shielding cutoff elevation angle in degrees.
    --height <height>  Shield height [default: 600].
    -o <figfile>      Output figure filename (PDF, PNG, etc.).
    -h --help         Show this help message.
"""
import matplotlib.pyplot as plt
import numpy as np
from lib_shield import bestrbt, trace_ray
from docopt import docopt


def main(T, h, outfile):
    # Single-stage shield geometry.
    # Design constraints: cutoff elevation T, shield height h, and a sufficiently wide base.
    # Near 60 degrees, the limiting ray grazes the wall and the base width approaches zero; no solution exists beyond that limit.

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    plt.subplots_adjust(bottom=0.2)

    if T < 59.5:
        rays_single(T, h, ax)
    else:
        ax.text(-150, 550, "No single case solution")

    # Two-stage shield geometry.
    if T > 53.5:
        rays_double(T, h, ax)

    ax.set_aspect('equal')
    ax.set_title(f"Shield single/double type for elevation angle {T} degree")
    if outfile:
        plt.savefig(outfile)

    plt.show()


# -----------------------------------------------------------------------------------------
def rays_double(T, h, ax):
    rapbest, rbtbest, r1best, h1best, Tdbest, raps, rbts = bestrbt(T, h)
    if ax is not None:
        model = {"nlines": 6, "start_id": 0, "aperture_id": 5, "bottom_id": 2,
                 0: {"A": np.array([-rapbest, h]), "B": np.array([-r1best, h1best])},
                 1: {"A": np.array([-r1best, h1best]), "B": np.array([-rbtbest, 0.])},
                 2: {"A": np.array([-rbtbest, 0.]), "B": np.array([rbtbest, 0])},
                 3: {"A": np.array([rbtbest, 0]), "B": np.array([r1best, h1best])},
                 4: {"A": np.array([r1best, h1best]), "B": np.array([rapbest, h])},
                 5: {"A": np.array([rapbest, h]), "B": np.array([-rapbest, h])}}
        x = [model[0]["A"][0]]
        y = [model[0]["A"][1]]
        for i in range(model["nlines"]):
            x.append(model[i]["B"][0])
            y.append(model[i]["B"][1])

        ax.plot(x, y, color='green')

        ax.text(0, -0.13,
                f"rap:{rapbest:.2f} rbt:{rbtbest:.0f} r1:{r1best:.2f} h1:{h1best:.2f} h:{h:.0f}",
                transform=ax.transAxes)
        for Tt in np.arange(T - 1.01, T + 1.01, 1):   # Trace rays near the limiting angle to illustrate reflections.
            P = model[0]["A"]
            U = np.array([np.cos(np.radians(Tt)), -np.sin(np.radians(Tt))])
            Ps, idy = trace_ray(model, P, U)
            if idy == model["aperture_id"]:
                ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='solid')
            else:
                ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='dotted')

        print(f"Td {Tdbest}")
        for Tt in np.arange(Tdbest - 1.01, Tdbest + 1.01, 2):   # Trace rays near the limiting angle to illustrate reflections.
            P = model[0]["A"]
            U = np.array([np.cos(np.radians(Tt)), -np.sin(np.radians(Tt))])
            Ps, idy = trace_ray(model, P, U)
            if idy == model["aperture_id"]:
                ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='solid')
            else:
                ax.plot(Ps[:, 0], Ps[:, 1], color='orange', linestyle='dotted')


# -----------------------------------------------------------------------------------------
def rays_single(T, h, ax):
    Tr = np.radians(T)   # Convert the cutoff angle to radians.
    # Geometric constraint: tan(T) = h / (r_ap + r_bt).
    rsum = h / np.tan(Tr)   # Sum of the aperture and base half-widths.
    print(f"rsum {rsum}")
    # Choose the wall slope so a ray striking its lower end reflects horizontally.
    # tan(T/2) = (r_ap - r_bt)/h
    rdif = np.tan(Tr / 2) * h
    print(f"rdif {rdif}")   # Diagnostic value for the original reference geometry.
    r_ap = (rdif + rsum) / 2
    print(f"r_ap {r_ap}")   # Diagnostic value for the original reference geometry.
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

        for Tt in np.arange(T - 1.01, T + 1.01, 1):   # Trace rays near the limiting angle to illustrate reflections.
            P = model[0]["A"]
            U = np.array([np.cos(np.radians(Tt)), -np.sin(np.radians(Tt))])
            Ps, idy = trace_ray(model, P, U)
            if idy == model["aperture_id"]:
                ax.plot(Ps[:, 0], Ps[:, 1], color='blue', linestyle='solid')
            else:
                ax.plot(Ps[:, 0], Ps[:, 1], color='blue', linestyle='dotted')

        ax.text(0, -0.17,
                f"rap:{r_ap:.2f} rbt:{r_bt:.0f} h:{h:.0f}",
                transform=ax.transAxes)


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = docopt(__doc__)
    T = float(args['--angle'])
    h = float(args['--height'])
    if not (0 < T < 90):
        raise SystemExit('Error: --angle must be between 0 and 90 degrees.')
    if h <= 0:
        raise SystemExit('Error: --height must be positive.')
    savefile = args['-o']
    main(T, h, savefile)
