# Spacecraft Attitude Calculation Example

## Overview

This example calculates the spacecraft attitude over one orbital revolution for six different
solar ecliptic longitudes.

The spacecraft follows a prescribed attitude profile consisting of an observation mode,
a non-observation mode, and transition intervals between them.
The observation target is the Galactic Center.

The calculation and visualization processes are separated:

1. `orb_att.py` calculates the orbit and attitude and saves the results in a compressed NumPy NPZ file.
2. `orb_view.py` reads the NPZ file and generates angular plots and three-dimensional attitude visualizations.
3. Additional numerical representations can be derived from the NPZ file without repeating the attitude calculation.

The NPZ file is the primary numerical output of this example.

For a description of the attitude model and its results,
see the [Spacecraft Attitude Modeling](https://indepex.github.io/projects/satellite/attitude/) page.
 *(Replace this URL with the published page URL.)*

## Directory Structure

```text
examples/attitude/
├── README.md
└── output/
    ├── attitude.npz
    ├── attitude_plot00.pdf
    ├── ...
    ├── attitude_plot05.pdf
    ├── attitude_view00.png
    ├── ...
    ├── attitude_view05.png
    ├── attitude_view00.html
    ├── ...
    └── attitude_view05.html
```

The input configuration files are stored in `examples/orbit/input/` and are shared with the orbital calculation example.

The output directory contains the calculated numerical data and visualization files.

## Simulation Conditions

| Parameter | Value |
|---|---|
| Orbital altitude | 570 km |
| Local time of ascending node (LTAN) | 06:00 |
| Observation target | Galactic Center |
| Target position-angle parameter | −2° |
| Solar ecliptic longitudes | 312°, 0°, 42°, 132°, 180°, 222° |
| Orbital revolutions | One per case |
| Sampling interval | 5° |

The six cases are stored together in a single NPZ file.

## Requirements

The example requires Python 3 and the following packages:

- NumPy
- Astropy
- sgp4
- docopt
- json5
- Matplotlib
- PyVista
- Cartopy

The required packages depend on the operation being performed. NumPy, Astropy, sgp4, docopt, and json5 are used in the calculation, while Matplotlib and PyVista are used for visualization.

Run the commands below from the repository root directory.

## Running the Calculation

Generate the attitude data for all six cases:

```bash
python orb_att.py \
    -s examples/attitude/output/attitude.npz \
    --orb examples/orbit/input/orbit.json5 \
    --obj examples/orbit/input/GC_pos.json5 \
    --sln 312,0,42,132.0,180,222.0 \
    --ltan 6 --alt 570 --tpa -2
```

The output is:

```text
examples/attitude/output/attitude.npz
```

The NPZ file contains the orbital and attitude data for all six cases.

The calculation is required only once. The visualization programs read the saved results.

## Generating Angular Plots

Generate the angular plots from the NPZ file:

```bash
python orb_view.py \
    examples/attitude/output/attitude.npz \
    -f \
    -d examples/attitude/output \
    -n attitude
```

Six PDF files are generated:

```text
attitude_plot00.pdf
attitude_plot01.pdf
attitude_plot02.pdf
attitude_plot03.pdf
attitude_plot04.pdf
attitude_plot05.pdf
```

Each plot shows the angular relationships between the spacecraft, the Sun, the observation target, and Earth as functions of orbital phase.

## Generating Three-Dimensional Views

Generate three-dimensional attitude visualizations:

```bash
python orb_view.py \
    examples/attitude/output/attitude.npz \
    -v -m arrow \
    -d examples/attitude/output \
    -n attitude
```

The arrow model represents the spacecraft body axes:

| Color | Spacecraft axis |
|---|---|
| Red | Z |
| Green | X |
| Blue | Y |

The spacecraft attitude is drawn at 5° orbital-phase intervals, with phase labels every 30°.

For each case, the program generates a PNG image and an interactive HTML file:

```text
attitude_view00.png
attitude_view00.html
...
attitude_view05.png
attitude_view05.html
```

The `-f` and `-v` options are mutually exclusive and must be used in separate commands.

## NPZ Data Format

The NPZ file contains the number of cases, the number of orbital-phase intervals, and the numerical results for each case.

Global entries include:

| Key | Description |
|---|---|
| `ncase` | Number of cases (6) |
| `npos` | Number of orbital-phase intervals (72) |
| `command_line` | Command-line arguments used for the calculation, if recorded |
| `coordinate_system` | Description of the coordinate-system approximations |

Each case uses a prefix such as `case000`, `case001`, ..., `case005`.

### Case Parameters

| Key suffix | Description |
|---|---|
| `_altitude` | Orbital altitude, km |
| `_ltan` | Local time of ascending node, hours |
| `_Torb` | Orbital period |
| `_target_name` | Observation target name |
| `_toObj` | Unit vector toward the observation target |
| `_year` | Simulation year |
| `_sun_longitude` | Solar ecliptic longitude, degrees |
| `_target_pa` | Target position-angle parameter, degrees |
| `_roll_angle` | Calculated roll correction, degrees |

### Time, Orbit, and Attitude

| Key suffix | Description | Expected shape |
|---|---|---|
| `_time_jd` | Sample times, Julian Date | `(73,)` |
| `_position` | Satellite position, TEME, km | `(73, 3)` |
| `_velocity` | Satellite velocity, TEME, km/s | `(73, 3)` |
| `_theta` | Orbital phase angles, degrees | `(73,)` |
| `_toSat` | Unit vector from Earth center toward the satellite | `(73, 3)` |
| `_toSun` | Unit vector toward the Sun | `(73, 3)` |
| `_SatX` | Spacecraft X-axis unit vector | `(73, 3)` |
| `_SatY` | Spacecraft Y-axis unit vector | `(73, 3)` |
| `_SatZ` | Spacecraft Z-axis unit vector | `(73, 3)` |
| `_jdc_jd` | Reference Julian Date | Scalar |

The 73 samples correspond to the 72 intervals covering one complete orbital revolution, including both endpoints.

### Reading the NPZ File

The data can be inspected using NumPy:

```python
import numpy as np

with np.load("examples/attitude/output/attitude.npz") as data:
    print("Number of cases:", data["ncase"].item())
    print("Number of intervals:", data["npos"].item())

    for k in range(int(data["ncase"])):
        p = f"case{k:03d}"
        print(
            k,
            "solar longitude =", data[f"{p}_sun_longitude"].item(),
            "position shape =", data[f"{p}_position"].shape,
        )
```

The spacecraft axes can be retrieved independently:

```python
import numpy as np

with np.load("examples/attitude/output/attitude.npz") as data:
    sat_x = data["case000_SatX"]
    sat_y = data["case000_SatY"]
    sat_z = data["case000_SatZ"]

    print("X-axis at first sample:", sat_x[0])
    print("Y-axis at first sample:", sat_y[0])
    print("Z-axis at first sample:", sat_z[0])
```

The three axis vectors define the spacecraft orientation at each sample.

They can be used to construct rotation matrices, quaternions, Euler angles, or attitudes relative to other reference frames.

No particular Euler-angle convention is imposed by the stored representation.

## Coordinate Systems and Limitations

The satellite position and velocity are calculated using SGP4 and expressed in the TEME frame.

The Sun direction is obtained using Astropy, while the target direction is specified in ICRS coordinates. These directions are treated as belonging to an approximately common inertial coordinate system.

The differences between these coordinate frames are neglected in this preliminary analysis.

The attitude profile is kinematic. It does not model spacecraft inertia, actuator performance, control errors, or the detailed dynamics of attitude maneuvers.

Transition attitudes are obtained by quaternion interpolation over prescribed orbital-phase intervals.

The numerical output is suitable for examining geometric relationships and for developing further analysis methods, subject to these approximations.
