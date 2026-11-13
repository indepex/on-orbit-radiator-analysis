# On-Orbit Spacecraft Radiator Analysis

## Overview

This repository provides Python and C tools for evaluating the radiative cooling performance of a spacecraft radiator used to cool payload detectors. The analysis accounts for spacecraft geometry, shielding, orbital motion, and attitude-dependent exposure to the Earth and Sun.

The workflow combines:

- Generation of synthetic Two-Line Element sets (TLEs) and orbit propagation using SGP4.
- Calculation of prescribed spacecraft attitudes and the directions of the Sun, Earth, and observation target.
- Evaluation of orbit-dependent Earth-radiation coupling using a simplified dice (cube-shaped) spacecraft model.
- Two-dimensional studies of thermal-radiation shield geometry.
- Three-dimensional ray tracing of incident radiation, reflection, and absorption.
- Calculation of directional radiation emitted from a radiator and its coupling to the Earth along an orbit.
- Preliminary radiative heat-budget estimates for detector cooling.

Numerical calculations and visualization are separated where practical. Intermediate results, including NPZ datasets, can be reused without repeating computationally intensive steps.

The software is intended for preliminary engineering studies rather than detailed spacecraft thermal verification.

## Examples

The examples include input data, representative outputs, and instructions for reproducing the calculations. **Run the documented commands from the repository root.**

| Example | Description |
| --- | --- |
| [Satellite orbit](examples/orbit/README.md) | Generate a synthetic TLE for a Sun-synchronous orbit and plot its ground track. |
| [Spacecraft attitude](examples/attitude/README.md) | Calculate orbit-dependent spacecraft attitude and visualize its geometry relative to the Sun, Earth, and observation target. |
| [Earth radiation: dice model](examples/thermal/dice/README.md) | Examine orbital-phase, seasonal, and attitude dependence of Earth coupling using a simplified cube-shaped satellite; the results inform shield design. |
| [Two-dimensional thermal shields](examples/thermal/shield-2d/README.md) | Compare shield configurations and examine the geometric trade-offs in suppressing unwanted radiation. |
| [Three-dimensional radiator analysis](examples/thermal/3d-model/README.md) | Model spacecraft shielding, trace incident sunlight, calculate directional radiator emission and Earth coupling, and estimate a radiative heat budget. |


## Calculation Workflow

The radiator design is developed in stages. A simplified dice (cube-shaped) satellite model is first used to examine how Earth-radiation coupling varies with orbital phase, season, and spacecraft attitude. This helps identify unfavorable viewing directions and motivates the subsequent shield design. The two-dimensional shield study then explores suitable geometries before the detailed three-dimensional analysis.

```text
Orbit configuration
       |
       +--> get_tle.py --> orbit.tle --> orbit_map.py --> ground-track map
       |
       +--> orb_att.py --> attitude.npz --> orb_view.py --> attitude plots/views
                                  |
                                  v
                       Simplified dice satellite
                       (earth_coupling.py)
                                  |
                                  v
                       Earth-coupling patterns
                       versus orbit and season
                                  |
                                  v
                       2D shield design study
                       (shield2d/ tools)
                                  |
                                  v
                       3D spacecraft geometry
                                  |
                 +----------------+----------------+
                 |                |                |
          model_view.py      from_far.py     outward_emission.py
          model image       Sun incidence         |
                                            outward.npz
                                                 |
                             attitude.npz -------+
                                                 |
                                             earth_sat.py
                                                 |
                                                 v
                                        Earth-coupling curves
                                                 |
                                                 v
                                        Preliminary heat budget
```

The dice-model coupling plots are included in `examples/thermal/dice/output/`. The 2D shield results are described in [the shield example](examples/thermal/shield-2d/README.md), and the detailed 3D calculations in [the radiator example](examples/thermal/3d-model/README.md). The arrows indicate the engineering design sequence; not every stage is a direct file-to-file dependency.

## Repository Structure

```text
on-orbit-radiator-analysis/
├── README.md
├── LICENSE                 # MPL-2.0: source code and configuration
├── LICENSE-DOCS            # CC0-1.0: documentation and original figures
├── Makefile
├── get_tle.py, orbit_map.py
├── orb_att.py, orb_view.py
├── earth_coupling.py, earth_sat.py
├── model_view.py, from_far.py
├── outward_emission.py, outward_view.py
├── lib_*.py
├── lib_radcouple.c, lib_radcouple.h
├── lib_vec.c, lib_vec.h
├── shield2d/
│   └── Python tools for two-dimensional shield analysis
└── examples/
    ├── orbit/
    ├── attitude/
    └── thermal/
        ├── dice/
        ├── shield-2d/
        └── 3d-model/
```

The source programs are located primarily in the repository root, with the two-dimensional shield tools under `shield2d/`. Each example contains its own documentation and relevant input and output files. Some output datasets are included to allow later analysis steps to be reproduced directly.

## Software Requirements and Build

The software uses Python 3 and a native C shared library for ray tracing. The following Python packages are imported by the supplied source files:

| Package | Purpose |
| --- | --- |
| NumPy | Numerical arrays and NPZ datasets |
| SciPy | Interpolation, filtering, and numerical optimization |
| Astropy | Astronomical coordinates and time calculations |
| sgp4 | Orbit propagation |
| numpy-quaternion | Spacecraft attitude rotations (`import quaternion`) |
| docopt | Command-line argument parsing |
| json5 | Configuration files |
| Matplotlib | Two-dimensional plotting |
| Cartopy | Ground-track maps |
| PyVista | Three-dimensional visualization |
| openpyxl | Spreadsheet output and processing |
| cycler | Plot style cycles |

The source also imports `erfa` (provided by the `pyerfa` package, normally installed as an Astropy dependency). To run the `Makefile`'s Python source checks, install `flake8` as well.

For example, the Python dependencies can be installed with:

```bash
python3 -m pip install numpy scipy astropy sgp4 numpy-quaternion docopt json5 \
    matplotlib cartopy pyvista openpyxl cycler flake8
```

The native library is built from `lib_radcouple.c` and `lib_vec.c` using GCC with OpenMP support. On an Ubuntu system, install the C build tools if needed:

```bash
sudo apt install build-essential
```

From the repository root, build the library with:

```bash
make lib_radcouple.so
```

This generates `lib_radcouple.so` alongside the Python sources. The Python wrapper `lib_radcpl.py` loads that file using `ctypes`.

Alternatively, running:

```bash
make
```

performs **both** the shared-library build and `flake8` checks on the Python files in the repository root. The `check` target does not currently include the Python files under `shield2d/`.

These instructions describe the supplied `Makefile`; they do not imply that the examples have been executed or tested on every platform.

## Running the Examples

Start with the [orbit example](examples/orbit/README.md) and [attitude example](examples/attitude/README.md) to understand the orbit-and-attitude dataset. The [dice satellite model example](examples/thermal/dice/README.md) evaluates the orbital and seasonal dependence of Earth-radiation coupling and provides a basis for radiator shielding design. The [2D shield example](examples/thermal/shield-2d/README.md) then examines shielding concepts, while the [3D radiator example](examples/thermal/3d-model/README.md) combines geometric ray tracing and orbit-dependent radiative coupling.

For example, the three-dimensional analysis reuses `examples/attitude/output/attitude.npz` together with `examples/thermal/3d-model/output/outward.npz`. Both datasets are supplied as representative outputs.

## Assumptions and Scope

The examples use prescribed orbit and attitude conditions and simplified geometric and radiative models. In particular:

- The orbit examples use synthetic TLEs and SGP4 propagation.
- The attitude profiles are kinematic prescriptions, not simulations of the attitude-control system.
- Ray-tracing results depend on the modeled geometry, surface behavior, and sampling resolution.
- Sun avoidance should be checked over the solar incidence angles expected during intended operations.
- Conversion of geometric Earth coupling into heat loads requires assumptions about infrared radiation and surface properties.
- The preliminary heat budgets do not include a complete transient thermal model, all conductive paths, or detailed temperature distributions.

The limitations and assumptions specific to each calculation are documented in its example README.

## License

The Python and C source code, header files, and configuration files in this repository are licensed under the **Mozilla Public License, Version 2.0 (MPL-2.0)**. See [LICENSE](LICENSE). When distributing modified MPL-covered files, you must make the corresponding source code available under the terms of the MPL-2.0. Separately developed files may be distributed under different terms, subject to the license conditions.

The original documentation, including this README and the example README files, original explanatory figures, and generated calculation results (including NPZ and CSV datasets, plots, and images) are dedicated to the public domain under **CC0 1.0 Universal**. See [LICENSE-DOCS](LICENSE-DOCS). These materials may be copied, modified, and redistributed, including commercially, without attribution. This dedication applies only to materials for which the project has the necessary rights; third-party materials, if any, retain their own terms.

The software is provided without warranty, as specified in the MPL-2.0.

## Project Website

This repository is part of the [indepex scientific and engineering projects](https://indepex.github.io/) collection. The website provides broader system-level context, while this repository contains the analysis programs, example data, and reproducible calculations.
