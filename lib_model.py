import numpy as np
import pyvista as pv
import json5


# Geometry type identifiers (consistent with lib_radcouple.h).
RECT, CYLIND, CUBE, TRI, DISC = 0, 1, 2, 3, 4

HIT_NONE = 0
HIT_ABSORB = 1
HIT_REFLECT = 2


# -----------------------------------------------------------------------------------------
def Arrow_Satellite():
    # Create a satellite model represented by its body-axis arrows.
    larrow = 800000
    C, X, Y, Z = np.array([0, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])
    Xarrow = pv.Arrow(start=C - (larrow / 2) * X, direction=X, scale=larrow)
    Yarrow = pv.Arrow(start=C - (larrow / 2) * Y, direction=Y, scale=larrow)
    Zarrow = pv.Arrow(start=C - (larrow / 2) * Z, direction=Z, scale=larrow)
    cmap = ['red', 'green', 'blue']
    opac = [1, 1, 1]
    mblock = pv.MultiBlock([Zarrow, Xarrow, Yarrow])

    return mblock, cmap, opac


# -----------------------------------------------------------------------------------------
def flagname(flg):
    if flg == 0:
        return "HIT_NONE"
    elif flg == 1:
        return "HIT_ABSORB"
    elif flg == 2:
        return "HIT_REFLECT"
    else:
        return f"Unknown flag {flg}"


# -----------------------------------------------------------------------------------------
def id2name(model, idx):
    # Return the name of the surface with the specified index.
    if idx >= 0:
        return model[idx]['name']
    return "none"


# -----------------------------------------------------------------------------------------
def name2id(model, name):
    # Return the index of the first surface matching the given name.
    for i, m in enumerate(model):
        if m["name"] == name:
            return i
    raise ValueError(f"surface name not found: {name}")


# -----------------------------------------------------------------------------------------
def generate_arrows(model):
    # Create surface-normal arrows for inspecting the model geometry.
    block = []
    cmap = []
    opac = []

    for item in model:
        if item["type"] == RECT:
            p0, p1, p2, p3 = item['p0'], item['p1'], item['p2'], item['p3']
            s = (p0 + p1 + p2 + p3) / 4
            d = np.cross((p1 - p0), (p3 - p0))
            # size = min(np.linalg.norm(p1 - p0), np.linalg.norm(p3 - p0))
            size = 300
            block.append(pv.Arrow(start=s, direction=d, scale=size))
            cmap.append('red')
            opac.append(1)
        elif item["type"] == TRI:
            p0, p1, p2 = item['p0'], item['p1'], item['p2']
            s = (p0 + p1 + p2) / 3
            d = np.cross((p1 - p0), (p2 - p0))
            # size = min(np.linalg.norm(p1 - p0), np.linalg.norm(p3 - p0))
            size = 300
            block.append(pv.Arrow(start=s, direction=d, scale=size))
            cmap.append('red')
            opac.append(1)

    return block, cmap, opac


# -----------------------------------------------------------------------------------------
def generate_block(model):
    """
    Return PyVista geometry, colors, and opacities for the model.

    Returns:
    -------
    tuple

        block : list of pyvista.DataSet
            Geometry objects representing the model surfaces.
        cmap : list
            Color assigned to each geometry object.
        opac : list
            Opacity assigned to each geometry object.
    """
    block = []
    cmap = []
    opac = []

    for item in model:
        if item["type"] == RECT:
            blk = rect(item)
        elif item["type"] == CYLIND:
            blk = cylinder(item)
        elif item["type"] == CUBE:
            blk = cube(item)
        elif item["type"] == TRI:
            blk = tri(item)
        elif item["type"] == DISC:
            blk = disc(item)

        else:
            raise ValueError(f"unknown type: {item['type']}")

        block.append(blk)

        # Render surfaces with a negative reflect flag as translucent.
        cmap.append(item["color"])
        if 'reflect' in item and item['reflect'] < 0:
            opac.append(0.1)
        else:
            opac.append(1.0)

    return block, cmap, opac


# -----------------------------------------------------------------------------------------
def cube(data):
    # Create a PyVista box from the specified dimensions and center.
    return pv.Cube(center=data['center'], x_length=data['x_len'],
                   y_length=data['y_len'], z_length=data['z_len'])


# -----------------------------------------------------------------------------------------
def cylinder(data):
    # Create an uncapped cylindrical surface.
    return pv.Cylinder(center=data['p0'], direction=data['p1'],
                       radius=data['radius'], height=data['height'], capping=False)


# -----------------------------------------------------------------------------------------
def disc(data):
    # Create a circular disk.
    return pv.Disc(center=data['p0'], normal=data['p1'], outer=data['radius'], c_res=50)


# -----------------------------------------------------------------------------------------
def rect(data):
    # Create a quadrilateral surface from four vertices.
    p0 = data['p0']
    p1 = data['p1']
    p2 = data['p2']
    p3 = data['p3']
    return pv.PolyData(np.array([p0, p1, p2, p3]), [4, 0, 1, 2, 3])


# -----------------------------------------------------------------------------------------
def tri(data):
    # Create a triangular surface from three vertices.
    p0 = data['p0']
    p1 = data['p1']
    p2 = data['p2']
    return pv.PolyData(np.array([p0, p1, p2]), [3, 0, 1, 2])


# -----------------------------------------------------------------------------------------
def load_payload_model(modelfile):
    """
    Load a payload model from JSON5 and construct surface and display geometry.
    Returns (data, cubes):
        data: ray-tracing surfaces, including a dummy starting surface.
        cubes: display-only geometry, excluded from ray tracing.
    """
    # Load the model geometry and parameters from the JSON5 file.
    with open(modelfile, 'r', encoding='utf-8') as f:
        rdata = json5.load(f)

    X, Y, Z = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])
    data = []   # Surfaces used in ray tracing.
    cubes = []  # Geometry used only for visualization.
    # Add a dummy surface for rays starting in free space; exclude it from physical intersections.
    p0 = 0.1 * X + 0.1 * Y + Z
    p1 = -0.1 * X + 0.1 * Y + Z
    p2 = -0.1 * X - 0.1 * Y + Z
    p3 = 0.1 * X - 0.1 * Y + Z
    d = setdata("Dummy", "darkgreen", 0, RECT, p0, p1, p2, p3, X, 0, 0)
    data.append(d)

    # Geometry excluded from ray tracing.
    # Use the top-center of the spacecraft bus as the coordinate origin.
    # Spacecraft bus.
    d = {"name": "Bus", "type": CUBE, "center": -rdata['bus_size'] / 2 * Z, "color": "gold",
         "x_len": rdata['bus_size'], "y_len": rdata['bus_size'], "z_len": rdata['bus_size']}
    cubes.append(d)
    # Solar array panels.
    d = {"name": "SapA", "type": CUBE, "center": rdata['bus_size'] / 2 * X + (rdata['SapW']
         + rdata['SapOff']) * Y - rdata['bus_size'] / 2 * Z, "color": "darkblue",
         "x_len": rdata['SapD'], "y_len": rdata['SapW'] * 2, "z_len": rdata['SapH']}
    cubes.append(d)
    d = {"name": "SapB", "type": CUBE, "center": rdata['bus_size'] / 2 * X - (rdata['SapW']
         + rdata['SapOff']) * Y - rdata['bus_size'] / 2 * Z, "color": "darkblue",
         "x_len": rdata['SapD'], "y_len": rdata['SapW'] * 2, "z_len": rdata['SapH']}
    cubes.append(d)
    # Lower the telescope housing slightly to avoid overlap with its absorbing top surface.
    d = {"name": "TelBox", "type": CUBE, "center": (rdata['tbox_height'] / 2 - 0.1) * Z,
         "color": "green",
         "x_len": rdata['tbox_width'], "y_len": rdata['tbox_width'],
         "z_len": rdata['tbox_height']}
    cubes.append(d)
    # Telescope housing top [A]: absorbing on both sides (reflect = 0).
    tw, te, th = rdata['tbox_width'], rdata['tbox_extend'], rdata['tbox_height']
    A0   = tw / 2       * X + tw / 2 * Y + th * Z
    A1 = -(tw / 2 + te) * X + tw / 2 * Y + th * Z
    A2 = -(tw / 2 + te) * X - tw / 2 * Y + th * Z
    A3   = tw / 2       * X - tw / 2 * Y + th * Z
    d = setdata("TelBoxTop", "darkgreen", 0, RECT, A0, A1, A2, A3, X, 0, 0)
    data.append(d)
    # Telescope baffle [B]: reflective on the outside, absorbing on the inside (reflect = 1).
    B0 = - rdata['baffle_offset'] * X + (th + rdata['baffle_height'] / 2) * Z
    d = setdata("TelBaffle", "#FFF4DB", 1, CYLIND, B0, Z, X, X, X,
                rdata['baffle_height'], rdata['baffle_dia'] / 2)
    data.append(d)
    # Absorbing disk inside the telescope baffle.
    B1 = - rdata['baffle_offset'] * X + (th + rdata['baffle_height'] * 0.8) * Z
    d = setdata("TelVane", "#606060", 0, DISC, B1, Z, X, X, X, 0, rdata['baffle_dia'] / 2)
    data.append(d)

    # Spacecraft envelope [C] (visualization geometry).
    C0 = rdata['env_h1'] / 2 * Z
    d = setdata("Envelope1", "#FF0000", -1, CYLIND, C0, Z, X, X, X,
                rdata['env_h1'], rdata['env_d1'] / 2)
    cubes.append(d)
    height = rdata['env_h3'] - rdata['env_h2']
    h, dh = 0, 50
    rt, rb = rdata['env_d4'] / 2, rdata['env_d3'] / 2
    while h <= height:
        r = (h + dh / 2) / height * (rt - rb) + rb
        C1 = (h + dh / 2 + rdata['env_h2']) * Z
        d = setdata("Envelope", "#FF0000", -1, CYLIND, C1, Z, X, X, X, dh, r)
        cubes.append(d)
        h += dh

    # Radiator [D]: absorbing on both sides (reflect = 0).
    w1, w2 = rdata['rad_width_1'], rdata['rad_width_2']
    d1, d2, h1 = rdata['rad_offset'], rdata['rad_depth'], rdata['rad_h1']
    D0 = -d1 * X + w1 / 2 * Y + (th + h1) * Z
    D1 = -(d1 + d2) * X + w2 / 2 * Y + (th + h1) * Z
    D2 = -(d1 + d2) * X - w2 / 2 * Y + (th + h1) * Z
    D3 = -d1 * X - w1 / 2 * Y + (th + h1) * Z
    data.append(setdata("Radiator", "#202020", 0, RECT, D0, D1, D2, D3, X, 0, 0))
    # First-stage radiator shield [E].
    w3, w4 = rdata['rad_width_3'], rdata['rad_width_4']
    d3, d4, h2 = rdata['rad_dx1'], rdata['rad_dx2'], rdata['rad_h2']
    E0 = -(d1 - d3) * X + w3 / 2 * Y + (th + h1 + h2) * Z
    E1 = -(d1 + d2 + d4) * X + w4 / 2 * Y + (th + h1 + h2) * Z
    E2 = -(d1 + d2 + d4) * X - w4 / 2 * Y + (th + h1 + h2) * Z
    E3 = -(d1 - d3) * X - w3 / 2 * Y + (th + h1 + h2) * Z
    data.append(setdata("RadRef1", "#FFF4DB", 3, RECT, E0, E1, D1, D0, X, 0, 0))
    data.append(setdata("RadRef2", "#FFF4DB", 3, RECT, E2, E3, D3, D2, X, 0, 0))
    data.append(setdata("RadRef3", "#FFF4DB", 3, RECT, E3, E0, D0, D3, X, 0, 0))
    data.append(setdata("RadRef4", "#FFF4DB", 3, RECT, E1, E2, D2, D1, X, 0, 0))
    # Second-stage radiator shield [F].
    w5, w6 = rdata['rad_width_5'], rdata['rad_width_6']
    d5, d6, h3 = rdata['rad_dx3'], rdata['rad_dx4'], rdata['rad_h3']
    F0 = -(d1 - d5) * X + w5 / 2 * Y + (th + h1 + h3) * Z
    F1 = -(d1 + d2 + d6) * X + w6 / 2 * Y + (th + h1 + h3) * Z
    F2 = -(d1 + d2 + d6) * X - w6 / 2 * Y + (th + h1 + h3) * Z
    F3 = -(d1 - d5) * X - w5 / 2 * Y + (th + h1 + h3) * Z
    data.append(setdata("RadRef5", "#FFF4DB", 3, RECT, F0, F1, E1, E0, X, 0, 0))
    data.append(setdata("RadRef6", "#FFF4DB", 3, RECT, F2, F3, E3, E2, X, 0, 0))
    data.append(setdata("RadRef7", "#FFF4DB", 3, RECT, F3, F0, E0, E3, X, 0, 0))
    # First-stage sunshield [G]–[H].
    sd1, sw1, sh1 = rdata['sld_d1'], rdata['sld_w1'], rdata['sld_h1']
    sd2, sw2 = rdata['sld_d2'], rdata['sld_w2']
    G0 = (tw / 2 + sd1) * X + sw1 / 2 * Y + th * Z
    G1 = (tw / 2 + sd1) * X - sw1 / 2 * Y + th * Z
    G2 = (tw / 2 - sd2) * X + sw2 / 2 * Y + th * Z
    G3 = (tw / 2 - sd2) * X - sw2 / 2 * Y + th * Z
    H0 = (tw / 2 + sd1) * X + sw1 / 2 * Y + (th + sh1) * Z
    H1 = (tw / 2 + sd1) * X - sw1 / 2 * Y + (th + sh1) * Z
    H2 = (tw / 2 - sd2) * X + sw2 / 2 * Y + (th + sh1) * Z
    H3 = (tw / 2 - sd2) * X - sw2 / 2 * Y + (th + sh1) * Z
    data.append(setdata("Sld1",  "#FFF4DB", 3, RECT, G0, G1, H1, H0, X, 0, 0))
    cubes.append(setdata("Sld1", "#FFFF00", 3, RECT, G0 + X, G1 + X, H1 + X, H0 + X, X, 0, 0))
    data.append(setdata("Sld2",  "#FFF4DB", 3, RECT, G2, G0, H0, H2, X, 0, 0))
    cubes.append(setdata("Sld2", "#FFFF00", 3, RECT, G2 + Y, G0 + Y, H0 + Y, H2 + Y, X, 0, 0))
    data.append(setdata("Sld3",  "#FFF4DB", 3, RECT, G1, G3, H3, H1, X, 0, 0))
    cubes.append(setdata("Sld3", "#FFFF00", 3, RECT, G1 - Y, G3 - Y, H3 - Y, H1 - Y, X, 0, 0))
    # Second-stage sunshield [I].
    sd3, sw3, sh2, sd4 = rdata['sld_d3'], rdata['sld_w3'], rdata['sld_h2'], rdata['sld_d4'],
    I0 = (tw / 2 - sd3) * X + sw3 / 2 * Y + (th + sh2) * Z
    I1 = (tw / 2 - sd3) * X - sw3 / 2 * Y + (th + sh2) * Z
    data.append(setdata("Sld4",  "#FFF4DB", 3, RECT, H0, H1, I1, I0, X, 0, 0))
    cubes.append(setdata("Sld4", "#FFFF00", 3, RECT, H0 + X, H1 + X, I1 + X, I0 + X, X, 0, 0))
    # I2 = I0 + k (H2-H0)  I2[0] = (tw / 2 - sd4) =  k (H2[0] - H0[0])
    I2 = I0 - (tw / 2 - sd4) / (H2[0] - H0[0]) * (H2 - H0)
    I3 = I1 - (tw / 2 - sd4) / (H3[0] - H1[0]) * (H3 - H1)
    # print('H21',H2[0]-H0[0],H2[0],H0[0])
    data.append(setdata("Sld5",  "#FFF4DB", 3, RECT, H2, H0, I0, I2, X, 0, 0))
    cubes.append(setdata("Sld5", "#FFFF00", 3, RECT, H2 + Y, H0 + Y, I0 + Y, I2 + Y, X, 0, 0))
    data.append(setdata("Sld6",  "#FFF4DB", 3, RECT, H1, H3, I3, I1, X, 0, 0))
    cubes.append(setdata("Sld6", "#FFFF00", 3, RECT, H1 - Y, H3 - Y, I3 - Y, I1 - Y, X, 0, 0))
    # Third-stage sunshield [J].
    sh3 = rdata['sld_h3']
    J0 = I0 + (sh3 - sh2) * Z
    J1 = I1 + (sh3 - sh2) * Z
    J2 = I2 + (sh3 - sh2) * Z
    J3 = I3 + (sh3 - sh2) * Z
    data.append(setdata("Sld7",  "#202020", 2, RECT, I0, I1, J1, J0, X, 0, 0))
    cubes.append(setdata("Sld7", "#FFFF00", 2, RECT, I0 + X, I1 + X, J1 + X, J0 + X, X, 0, 0))
    data.append(setdata("Sld8",  "#202020", 2, RECT, I2, I0, J0, J2, X, 0, 0))
    cubes.append(setdata("Sld8", "#FFFF00", 2, RECT, I2 + Y, I0 + Y, J0 + Y, J2 + Y, X, 0, 0))
    data.append(setdata("Sld9",  "#202020", 2, RECT, I1, I3, J3, J1, X, 0, 0))
    cubes.append(setdata("Sld9", "#FFFF00", 2, RECT, I1 - Y, I3 - Y, J3 - Y, J1 - Y, X, 0, 0))
    data.append(setdata("Slda",  "#202020", 2, TRI,  H3, J3, I3,  X, X, 0, 0))
    cubes.append(setdata("Slda", "#FFFF00", 2, TRI,  H3 + X, J3 + X, I3 + X, X, X, 0, 0))
    data.append(setdata("Sldb",  "#202020", 2, TRI,  I2, J2, H2, X, X, 0, 0))
    cubes.append(setdata("Sldb", "#FFFF00", 2, TRI,  I2 + X, J2 + X, H2 + X, X, X, 0, 0))

    return data, cubes


# -----------------------------------------------------------------------------------------
def setdata(name, color, alp, tp, p0, p1, p2, p3, p4, h, r):
    return {"name": name, "color": color, "reflect": alp, "type": tp,
            "p0": p0, "p1": p1, "p2": p2, "p3": p3, "p4": p4, "height": h, "radius": r}
