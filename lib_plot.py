import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from openpyxl import Workbook, load_workbook
import re


# -----------------------------------------------------------------------------------------
def plot_angle(data, plot_file):
    # Extract orbit and attitude data.
    toSun, toSat = data["toSun"], data["toSat"]
    altitude, ltan, roll = data["altitude"], data["ltan"], data["roll_angle"]
    # jdiso = data["jdiso"]
    slon = data["sun_longitude"]
    theta = data["theta"]
    toObj = data["toObj"]
    SatX, SatY, SatZ = data["SatX"], data["SatY"], data["SatZ"]
    # npos = toSat.shape[0] - 1
    # Set the plot title.
    title = f"Alt {altitude:.0f}km ltan {ltan:.2f} Roll {roll:.2f} S.lon {slon}"
    # Target elevation above the local horizontal plane.
    TargetElevations = np.degrees(np.arcsin(np.clip(toSat @ toObj, -1., 1.)))
    # Sun direction in the spacecraft body frame.
    sz = np.clip(np.sum(toSun * SatZ, axis=1), -1, 1)
    sx = np.clip(np.sum(toSun * SatX, axis=1), -1, 1)
    sy = np.clip(np.sum(toSun * SatY, axis=1), -1, 1)
    th = 90 - np.degrees(np.arccos(sz))
    ph = np.degrees(np.arctan2(sy, sx))

    # Set up the figure and axes.
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(2, 1, 1)
    ax3 = fig.add_subplot(2, 1, 2)
    # ax2 = ax3.twinx()

    ax1.set_ylim(-90, 90)
    ax1.plot(theta, TargetElevations, label="Target Elevation", color='gray', linestyle='-')
    ax1.plot(theta, th, label="Sun Elevation angle", color='magenta', linestyle='-')
    ax1.plot(theta, ph, label="Sun Azimuth angle", color='blue', linestyle='-')
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.85), ncol=5, fontsize=8, frameon=False)
    ax1.set_ylabel('Solar el./az. angle [degree]')

    ax3.set_xlabel('Satellite phase angle [degree]')
    ax3.set_ylim(-180, 180)

    cosX = np.sum(SatX * (-toSat), axis=1)
    cosY = np.sum(SatY * (-toSat), axis=1)
    cosZ = np.sum(SatZ * (-toSat), axis=1)

    EarthX = 90 - np.degrees(np.arccos(np.clip(cosX, -1.0, 1.0)))
    EarthY = 90 - np.degrees(np.arccos(np.clip(cosY, -1.0, 1.0)))
    EarthZ = 90 - np.degrees(np.arccos(np.clip(cosZ, -1.0, 1.0)))

    ax3.plot(theta, EarthX, label="Earth El. to X", color='magenta', linestyle='-')
    ax3.plot(theta, EarthY, label="Earth El. to Y", color='orange', linestyle='-')
    ax3.plot(theta, EarthZ, label="Earth El. to Z", color='red', linestyle='-')
    ax3.legend(loc='lower center', bbox_to_anchor=(0.5, 0.85), ncol=5, fontsize=8, frameon=False)
    ax3.set_ylabel('Earth elevation angle [degree]')

    ax1.set_title(title)
    if plot_file:
        plt.savefig(plot_file, bbox_inches="tight")
        plt.show()
    else:
        plt.show()


# -----------------------------------------------------------------------------------------
def plot_points(points, colors, modelfile, tot, absorbed, savefile, elevation=None, azimuth=None):
    """
    Plot points as a 3D scatter plot using the specified colors.
    """
    # Extract the Cartesian coordinates.
    x, y, z = points[:, 0],  points[:, 1],  points[:, 2]
    # Create the 3D scatter plot.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.scatter(x, y, z, c=colors, marker='o')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)

    modelname = Path(modelfile).name
    rabs = absorbed / tot * 100
    ax.set_title(
        f"Coupling with the Earth\n{modelname} total {tot:.0f}, abs {absorbed:.0f}({rabs:.1f}%)")

    if elevation is not None and azimuth is not None:
        ax.view_init(elev=elevation, azim=azimuth)

    def on_close(event):
        # Save the figure when its window closes.
        if savefile:
            fig.savefig(savefile)
            print(f"Saved figure to {savefile}.")

    # Register the window-close callback.
    fig.canvas.mpl_connect('close_event', on_close)

    plt.show()


# -----------------------------------------------------------------------------------------
def output_xlsx(data, xls_file, newfile=False):
    # Export orbit and spacecraft attitude data to an Excel workbook.
    # Extract the input data.
    theta = data["theta"]
    toSat = data["toSat"]
    toSun = data["toSun"]
    SatX, SatY, SatZ = data["SatX"], data["SatY"], data["SatZ"]
    # Sun direction in the spacecraft body frame.
    sx = np.clip(np.sum(toSun * SatX, axis=1), -1, 1)
    sy = np.clip(np.sum(toSun * SatY, axis=1), -1, 1)
    sz = np.clip(np.sum(toSun * SatZ, axis=1), -1, 1)
    # Earth direction in the spacecraft body frame.
    ex = np.clip(np.sum(-toSat * SatX, axis=1), -1, 1)
    ey = np.clip(np.sum(-toSat * SatY, axis=1), -1, 1)
    ez = np.clip(np.sum(-toSat * SatZ, axis=1), -1, 1)
    ezt = np.degrees(np.arccos(ez))
    # Derive the worksheet name from the Sun longitude.
    sheetname = "{}".format(data["sun_longitude"])
    # Replace characters prohibited in Excel worksheet names.
    sheetname = re.sub(r'[:\\/*?\[\]]', '_', sheetname)
    # Strip leading and trailing whitespace.
    sheetname = sheetname.strip()
    # Remove control characters, including newlines and tabs.
    sheetname = ''.join(ch for ch in sheetname if ord(ch) >= 32)
    sheetname = sheetname.strip().replace(' ', '_')
    sheetname = sheetname.replace('.', '_')
    # Create a workbook or append a new worksheet.
    if newfile:
        wb = Workbook()
        ws = wb.active
        ws.title = sheetname
    else:
        wb = load_workbook(xls_file)
        ws = wb.create_sheet(title=sheetname)

    # Set column headings.
    ws["A1"] = "Angle from Start [deg]"
    ws["B1"] = "Sat Position Unit Vector"
    ws["B2"], ws["C2"], ws["D2"] = "X", "Y", "Z"
    ws["E1"] = "Sun direction in Sat"
    ws["E2"], ws["F2"], ws["G2"] = "X", "Y", "Z"
    ws["H1"] = "Earth direction in Sat"
    ws["H2"], ws["I2"], ws["J2"] = "X", "Y", "Z"
    ws["K1"] = "Sat X-Axis in Ref"
    ws["K2"], ws["L2"], ws["M2"] = "X", "Y", "Z"
    ws["N1"] = "Sat Y-Axis in Ref"
    ws["N2"], ws["O2"], ws["P2"] = "X", "Y", "Z"
    ws["Q1"] = "Sat Z-Axis in Ref"
    ws["Q2"], ws["R2"], ws["S2"] = "X", "Y", "Z"
    ws.merge_cells("B1:D1")
    ws.merge_cells("E1:G1")
    ws.merge_cells("H1:J1")
    ws.merge_cells("K1:M1")
    ws.merge_cells("N1:P1")
    ws.merge_cells("Q1:S1")
    for i in range(len(theta)):
        k = i + 3
        ws[f"A{k}"] = theta[i]
        ws[f"B{k}"] = toSat[i][0]
        ws[f"C{k}"] = toSat[i][1]
        ws[f"D{k}"] = toSat[i][2]

        ws[f"E{k}"] = sx[i]
        ws[f"F{k}"] = sy[i]
        ws[f"G{k}"] = sz[i]
        ws[f"H{k}"] = ex[i]
        ws[f"I{k}"] = ey[i]
        ws[f"J{k}"] = ez[i]

        ws[f"K{k}"] = SatX[i][0]
        ws[f"L{k}"] = SatX[i][1]
        ws[f"M{k}"] = SatX[i][2]
        ws[f"N{k}"] = SatY[i][0]
        ws[f"O{k}"] = SatY[i][1]
        ws[f"P{k}"] = SatY[i][2]
        ws[f"Q{k}"] = SatZ[i][0]
        ws[f"R{k}"] = SatZ[i][1]
        ws[f"S{k}"] = SatZ[i][2]

        ws[f"T{k}"] = ezt[i]

    wb.save(xls_file)

    return
