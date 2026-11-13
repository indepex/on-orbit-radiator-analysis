import numpy as np
import pyvista as pv
import warnings
from pathlib import Path
from astropy import constants as AC
import math
import quaternion
from sgp4.earth_gravity import wgs84

from lib_utils import normalize
from lib_model import Arrow_Satellite


# -------------------------------------------------------------------------------------
def view_sat(data, model, outfile, htmlfile):
    """
    Visualize the Earth, spacecraft orbit, and satellite in 3D.

    Parameters
    ----------
    data : dict
        Orbit and attitude calculation results.
    model : str
        Satellite geometry model (currently "arrow").
    outfile, htmlfile : str or False
        Output paths for an image and an interactive HTML scene, respectively.
        If neither is specified, display the scene interactively.
    """
    ORG = np.array([0, 0, 0])
    toSun, toSat = data["toSun"], data["toSat"]
    altitude, ltan, roll = data["altitude"], data["ltan"], data["roll_angle"]
    # jdiso = data["jdciso"]
    slon = data["sun_longitude"]
    SatX, SatZ = data["SatX"], data["SatZ"]
    theta = data["theta"]
    npos = toSat.shape[0] - 1
    Sun = toSun[npos // 4]  # Representative Sun direction near the midpoint of the observation period.
    # Set up the Sun as a light source and create the Earth model.
    sun = pv.Light(light_type='scene light',
                   position=AC.au.value * Sun,   # Position in meters.
                   focal_point=ORG, color='white')
    # Earth.
    Re = wgs84.radiusearthkm  # Earth radius [km].
    globe = pv.Sphere(radius=Re * 1000, center=ORG)  # Coordinates in meters.
    # Normal vector to the orbital plane.
    npos1, npos2, npos3 = 0, int(npos // 3), int(npos * 2 // 3)  # Select three points along the orbit.
    v1, v2 = toSat[npos2] - toSat[npos1], toSat[npos3] - toSat[npos1]
    v = normalize(np.cross(v2, v1))    # Normal vector to the orbital plane.
    # Satellite orbit.
    Sat = toSat * (Re + altitude) * 1000   # Coordinates in meters.
    SatL = toSat * (Re + altitude * 1.5) * 1000   # Positions for phase labels.
    orbit = pv.PolyData(Sat)
    orbit.lines = np.hstack([[len(Sat)], np.arange(len(Sat))])

    # Satellite geometry.
    if model == "arrow":
        SAT, SATc, SATo = Arrow_Satellite()
    else:
        raise ValueError(f"unknown model: {model}")

    """
    SAT, SATc, SATo = Satellite_block(model)
    """
    # Configure and render the scene.
    if outfile:  # Render off-screen when saving an image.
        plotter = pv.Plotter(off_screen=True)
    else:
        plotter = pv.Plotter()
    plotter.enable_parallel_projection()  # Use parallel projection.
    plotter.remove_all_lights()
    plotter.add_light(sun)
    plotter.add_mesh(globe, smooth_shading=True, ambient=0.0, diffuse=1.0,
                     specular=0.0, color=(0.20, 0.70, 1.00))
    plotter.add_mesh(orbit, color=(0.20, 0.70, 1.00), line_width=1, point_size=1)
    """
    # Arrow indicating the Sun direction.
    toSunArrow = pv.Arrow(start=sunlight * AC.R_earth.value * scale,
                          direction=sunlight, scale=5000)
    plotter.add_mesh(toSunArrow, color='red')
    """
    for i in range(0, npos, 1):
        # Render the satellite at the current orbital position.
        J = rotate_move_block(SAT, SatX[i], SatZ[i], Sat[i])
        for k in range(J.n_blocks):
            plotter.add_mesh(J[k], color=SATc[k], opacity=SATo[k])
        # Label orbital phase at 30-degree intervals.
        if i % 6 == 0:
            label_pos = np.asarray(SatL[i]) * 1.05
            plotter.add_point_labels(np.array([label_pos]), [f"{theta[i]:.0f}"],
                                     point_size=1, font_size=12, always_visible=True)
    """
    # Orbital axis.
    OrbAx = normalize(np.cross(Sat[int(npos / 2)], Sat[int(npos / 4)]))
    # Direction of the orbital rotation axis.
    OrbAxArrow = pv.Arrow(start=OrbAx * AC.R_earth.value * scale,
                          direction=OrbAx, scale=5000)
    plotter.add_mesh(OrbAxArrow, color='green')
    """
    # Fixed-position annotation.
    title = f"Alt {altitude:.0f}km\nltan {ltan:.2f}\nRoll {roll:.2f}\nSun lon. {slon}"
    plotter.add_text(title, position="upper_left", font_size=14)
    if model == "arrow":
        plotter.add_legend(
            labels=[["Satellite Z", "red"], ["Satellite X", "green"], ["Satellite Y", "blue"], ],
            bcolor="white", border=False, size=(0.18, 0.12), loc="lower right",)
    if ltan < 12:
        cpos = (v[0], v[1], v[2])
    else:
        cpos = (-v[0], -v[1], -v[2])
    plotter.camera_position = [cpos, (0, 0, 0), (0, 1, 0)]

    plotter.reset_camera()  # Reset the camera to fit the entire scene.
    plotter.camera.zoom(1.5)  # Zoom in (values above 1 enlarge; values below 1 reduce).

    if outfile:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            plotter.show(auto_close=False)
            plotter.screenshot(outfile)

    if htmlfile:
        plotter.export_html(htmlfile)

    if outfile or htmlfile:
        plotter.close()
    else:
        plotter.show()
        plotter.close()


# -----------------------------------------------------------------------------------------
def rotate_move_block(block, U, W, P):
    """
    Rotate and translate a PyVista MultiBlock into a new coordinate frame.

    Rotate the geometry about the origin so that its local X and Z axes align
    with U and W, respectively. Then translate the rotated geometry by P.

    Parameters
    ----------
    block : pyvista.MultiBlock
        Original 3D geometry.
    U, W : ndarray of shape (3,)
        Target directions of the local X and Z axes, respectively.
    P : ndarray of shape (3,)
        Translation vector specifying the new origin.

    Returns
    -------
    r_block : pyvista.MultiBlock
        Rotated and translated geometry.
    """
    NX = np.array([1, 0, 0])  # Original X axis.
    NZ = np.array([0, 0, 1])  # Original Z axis.

    # Allocate the transformed output blocks.
    r_block = pv.MultiBlock([None] * block.n_blocks)

    # Rotate the model so its local X and Z axes align with U and W.
    # The Y axis follows from the right-handed coordinate system.
    # First, rotate the Z axis onto W.
    rot = np.cross(NZ, W)
    rotlen = np.linalg.norm(rot)
    if not np.allclose(NZ, W, atol=1e-12):  # Skip the rotation if the axes already coincide.
        if rotlen > np.finfo(float).eps:
            rotax = rot / rotlen
        else:  # Use the X axis for a 180-degree rotation.
            rotax = NX
        theta = math.acos(np.clip(np.dot(NZ, W), -1.0, 1.0))
        # Apply the first rotation to the X axis.
        rot_quat = quaternion.from_rotation_vector(theta * rotax)
        NXr = quaternion.as_rotation_matrix(rot_quat).dot(NX)
    else:  # Axes already coincide; no rotation is needed.
        NXr = NX
        rot_quat = np.quaternion(1, 0, 0, 0)

    # Then rotate the transformed X axis about W by omega to align it with U.
    if not np.allclose(NXr, U, atol=1e-12):  # Skip the rotation if the axes already coincide.
        omega = math.atan2(np.dot(W, np.cross(NXr, U)), np.dot(NXr, U))  # rad, signed
        rot2_quat = quaternion.from_rotation_vector(omega * W)
    else:
        rot2_quat = np.quaternion(1, 0, 0, 0)
    # Compose the two quaternions to obtain the complete rotation.
    rotall_quat = rot2_quat * rot_quat
    # Convert the combined quaternion to an axis-angle rotation.
    rotation_axis = quaternion.as_rotation_vector(rotall_quat)
    angle_rad = (np.linalg.norm(rotation_axis))
    rotation = rotation_axis / angle_rad if angle_rad > 0 else np.array([1, 0, 0])
    angle = math.degrees(angle_rad)
    for i in range(block.n_blocks):
        mesh = block[i].rotate_vector(rotation, angle, inplace=False)
        r_block[i] = mesh.translate(P, inplace=False)

    return r_block


# -------------------------------------------------------------------------------------
def view_block(block, cmap, opacity, savefile, comment=None, campos=None, viewan=None, pscale=None):
    """
    Render a MultiBlock with per-block colors and opacities.

    Save the scene as an image or HTML file when a path is supplied;
    otherwise, display it interactively.

    Parameters
    ----------
    block : pyvista.MultiBlock
        Geometry to render.
    cmap : list
        Color assigned to each block.
    opacity : list
        Opacity assigned to each block.
    savefile : str or None
        Output image or HTML path; None displays the scene interactively.
    comment : str or None
        Text annotation fixed to the upper-left corner.
    campos : optional
        Camera position passed to PyVista.
    viewan : float or None
        Camera view angle.
    pscale : float or None
        Parallel projection scale.
    """
    # Configure the plotter.
    plotter = pv.Plotter()
    plotter.enable_parallel_projection()

    # Render each mesh block.
    for i in range(block.n_blocks):
        plotter.add_mesh(block[i], color=cmap[i], opacity=opacity[i],
                         render_points_as_spheres=False, point_size=0)
    if comment is not None:
        plotter.add_text(comment, position="upper_left", font_size=10)

    if campos is not None:
        plotter.camera_position = campos
    if viewan is not None:
        plotter.camera.view_angle = viewan
    if pscale is not None:
        plotter.camera.parallel_scale = pscale

    plotter.camera_set = True

    # Save the visualization when an output file is specified.
    if savefile:
        if Path(savefile).suffix == '.html':
            plotter.export_html(savefile)
            plotter.close()
        else:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                plotter.show(auto_close=False)
                plotter.screenshot(savefile)
                plotter.close()
    else:
        plotter.show()
        saved_cpos = plotter.camera_position
        saved_view_angle = plotter.camera.view_angle
        saved_scale = plotter.camera.parallel_scale
        print("cops ", saved_cpos)
        print("view angle", saved_view_angle)
        print("parallel scale", saved_scale)


# -----------------------------------------------------------------------------------------
def set_cam(azimuth, elevation, distance, fp_x, fp_y, fp_z):
    """
    Compute a PyVista camera position from azimuth and elevation.

    Parameters
    ----------
    azimuth, elevation : float
        Camera direction angles in degrees.
    distance : float
        Distance from the focal point to the camera offset.
    fp_x, fp_y, fp_z : float
        Focal point coordinates.

    Returns
    -------
    list
        Camera position, focal point, and view-up vector (positive Z).
    """
    azim = np.radians(azimuth)
    elev = np.radians(elevation)

    cx = distance * np.cos(elev) * np.cos(azim)
    cy = distance * np.cos(elev) * np.sin(azim)
    cz = distance * np.sin(elev)

    camera_position = [(cx, cy, cz), (fp_x, fp_y, fp_z), (0.0, 0.0, 1.0)]

    return camera_position
