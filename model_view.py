#!/usr/bin/env python3
r"""
Usage:
  model_view.py <model.json5> [-c] [-o <plotfile>]
  model_view.py (-h | --help)

Options:
  -h --help       Show this help.
  -c              Use the predefined camera.
  -o <plotfile>   Save the rendered plot to a file.
"""
import pyvista as pv
from docopt import docopt
from lib_view import set_cam, view_block
from lib_model import load_payload_model, generate_block, generate_arrows


def main(modelfile, savefile, cam, van, psc):

    # Load the payload model.
    model, cubes = load_payload_model(modelfile,)

    # Generate meshes for the bounding cubes, payload components, and axis arrows.
    cblk, cc, co = generate_block(cubes)
    block, cmap, opac = generate_block(model)
    arblk, arc, aro = generate_arrows(model)

    # Axis arrows are included in the combined mesh.

    # Combine all meshes into a single MultiBlock.
    combb = pv.MultiBlock(cblk + block + arblk)

    # Render the model or save the visualization.
    view_block(combb, cc + cmap + arc, co + opac + aro,
               savefile, campos=cam, viewan=van, pscale=psc)


if __name__ == '__main__':
    args = docopt(__doc__)
    modelfile = args['<model.json5>']
    savefile = args['-o']
    if args['-c']:
        cam = set_cam(125, 35, 18000, 0, 0, 150)  # Azimuth, elevation, distance, and focal-point coordinates (x, y, z).
        van = 40.                        # Camera view angle [deg].
        psc = 2000.                      # Parallel-projection scale.
    else:
        cam = None
        van = None
        psc = None

    main(modelfile, savefile, cam, van, psc)
