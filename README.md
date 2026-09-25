# Mission Analysis Project — Earth to Asteroid 3908 Nyx

Design and analysis of an interplanetary mission for the orbital insertion of a spacecraft around the near-Earth asteroid **3908 Nyx (1980 PA)**, starting from a Geostationary Transfer Orbit (GTO). The mission is designed within the **patched-conics approximation**, with every maneuver optimized to minimize the total propellant cost (Δv).

**Authors (Group 26):** [Francesco Bernardi](https://github.com/Bernaaaa), [Tancredi Bertolotti](https://github.com/ATOMIXXPOWER), [Giulia Anzivino](https://github.com/Giuly1111)

📄 Full report: [`docs/report/Report_Space_LAB_PoliMi.pdf`](docs/report/Report_Space_LAB_PoliMi.pdf)

---

## Overview

The mission is broken into three scenarios, each modeled as an independent two-body problem and later linked together via the Patched Conics method:

| # | Scenario | Description | Code |
|---|----------|-------------|------|
| 1 | **Reaching Parking Orbit** | Geocentric transfer from the initial GTO to the designated parking orbit; three alternative strategies are designed and compared | [`src/Project/Scenery 1`](src/Project/Scenery%201) |
| 2 | **Interplanetary Direct Transfer** | Heliocentric transfer from Earth's orbit to asteroid Nyx's orbit, solved via grid search + `fmincon` optimization | [`src/Project/Scenery 2`](src/Project/Scenery%202) |
| 3 | **Escape & Arrival Phases** | Hyperbolic escape from the Earth parking orbit and hyperbolic capture at Nyx, both coplanar and non-coplanar cases | [`src/Project/Scenery 3`](src/Project/Scenery%203) |

### Result at a glance

| Phase | Δv [km/s] | ΔT [10⁴ s] |
|---|---|---|
| Scenario 1 (bi-elliptic, optimal strategy) | 3.2254 | 206.188 |
| Scenario 2 (heliocentric transfer)* | 0** | 529.918 |
| Scenario 3 — departure (non-coplanar, optimal injection) | 3.8304 | 2.683 |
| Scenario 3 — arrival (Nyx capture) | 4.0171 | – |
| **Total** | **11.0729** | **738.789** |

\* Total time of flight: **85 days, 12 h, 11 min, 30 s**
\*\* The Scenario 2 Δv is already embedded in the hyperbolic excess velocities (v∞) used in Scenario 3; it is not summed separately to avoid double counting.

Optimization highlights:
- **Scenario 1:** comparing a standard bitangent strategy (7.53 km/s), a time-optimal strategy (8.88 km/s) and a bi-elliptic strategy found via a 500×500 grid search (3.23 km/s) → **~57% Δv saving** by selecting the bi-elliptic strategy.
- **Scenario 3:** comparing a fixed pericenter-injection departure (15.09 km/s) against a grid-search over the optimal injection point on the parking orbit (3.83 km/s, at θ ≈ 165°) → **~74% Δv saving**.

---

## Repository Structure

```
.
├── README.md
├── license
├── docs/
│   ├── data/            # Assignment sheets & numerical input data
│   ├── report/          # Final written report (PDF)
│   └── video/           # Trajectory animations (MP4)
└── src/
    ├── lab/             # Reusable astrodynamics toolbox (developed in the weekly labs)
    │   ├── L1/          # Coordinate conversions
    │   ├── L2/          # Orbital maneuvers
    │   └── L3/          # Time of flight
    └── Project/         # Mission-specific scripts, organized by scenario
        ├── Scenery 1/
        ├── Scenery 2/
        └── Scenery 3/
```

---

## `src/lab` — Core Astrodynamics Toolbox

Generic, reusable MATLAB functions developed across the lab sessions and used throughout the whole `Project` folder as the mission's computational backbone.

### L1 — Coordinate Conversions
| File | Description |
|---|---|
| [`par2car.m`](<src/lab/L1/par2car.m>) | Converts classical Keplerian elements (a, e, i, Ω, ω, θ) into Cartesian state vectors (**r**, **v**) |
| [`car2par.m`](<src/lab/L1/car2par.m>) | Inverse conversion: from Cartesian state vectors to Keplerian elements |

### L2 — Orbital Maneuvers
| File | Description |
|---|---|
| [`bitangentTransfer.m`](<src/lab/L2/bitangentTransfer.m>) | Computes a single Hohmann-like bitangent transfer between two coaxial orbits (Pericenter-Apocenter / Apocenter-Pericenter / etc.) |
| [`biellipticTransfer.m`](<src/lab/L2/biellipticTransfer.m>) | Computes a bi-elliptic transfer through an intermediate orbit, used to cut plane-change costs |
| [`changeOrbitalPlane.m`](<src/lab/L2/changeOrbitalPlane.m>) | Combined inclination/RAAN change maneuver at the optimal node, with full output/printing |
| [`changeOrbitalPlaneNoPrint.m`](<src/lab/L2/changeOrbitalPlaneNoPrint.m>) | Same as above, silent version (used inside optimization/grid-search loops) |
| [`changeOrbitalPlaneDeltaT.m`](<src/lab/L2/changeOrbitalPlaneDeltaT.m>) | Plane-change maneuver selecting the node reachable in minimum time, rather than minimum Δv |
| [`changePericenterArg.m`](<src/lab/L2/changePericenterArg.m>) | Argument-of-pericenter rotation maneuver (apsidal line rotation) at fixed a, e |

### L3 — Time of Flight
| File | Description |
|---|---|
| [`TOF_E.m`](<src/lab/L3/TOF_E.m>) | Time of flight via eccentric anomaly (elliptical orbits), with console output |
| [`TOF_E_NoPrint.m`](<src/lab/L3/TOF_E_NoPrint.m>) | Silent version, for use in loops/optimizations |
| [`TOF_M.m`](<src/lab/L3/TOF_M.m>) | Time of flight via mean anomaly / Kepler's equation, with console output |
| [`TOF_M_NoPrint.m`](<src/lab/L3/TOF_M_NoPrint.m>) | Silent version, for use in loops/optimizations |

---

## `src/Project` — Mission Scenarios

### Scenario 1 — Reaching Parking Orbit
`src/Project/Scenery 1/`

Transfers the spacecraft from the initial GTO (a = 24 400 km, e = 0.7283) to the final geocentric parking orbit (a = 42 145.32 km, e = 0.7659), requiring a combined plane change, apsidal rotation and shape change. Three strategies are compared: **Standard**, **ΔT-optimal** and **Δv-optimal (bi-elliptic)**.

| File | Description |
|---|---|
| [`scenery1_standard.m`](<src/Project/Scenery 1/scenery1_standard.m>) | Builds Strategy 1 & 2 (standard and time-optimal bitangent sequences) and computes their Δv/ΔT |
| [`scenery1_deltaV.m`](<src/Project/Scenery 1/scenery1_deltaV.m>) | Builds Strategy 3: bi-elliptic transfer through an optimized intermediate orbit |
| [`scenery1_deltaT.m`](<src/Project/Scenery 1/scenery1_deltaT.m>) | Time-of-flight computation/comparison across the three strategies |
| [`scenery1_optimization_plot.m`](<src/Project/Scenery 1/scenery1_optimization_plot.m>) | Grid search over the intermediate orbit's (rₐ, e) for the bi-elliptic strategy + Δv contribution plots (Fig. 1.4–1.5 in the report) |
| [`scenery1_plot_standard.m`](<src/Project/Scenery 1/scenery1_plot_standard.m>) | 3D trajectory plot for Strategies 1 & 2 |
| [`scenery1_plot_deltaV.m`](<src/Project/Scenery 1/scenery1_plot_deltaV.m>) | 3D trajectory plot for Strategy 3 (bi-elliptic) |
| [`istogram_comparison_s1.m`](<src/Project/Scenery 1/istogram_comparison_s1.m>) | Bar-chart comparison of Δv/ΔT across the three strategies (Fig. 1.6) |
| [`earth_texture.png`](<src/Project/Scenery 1/earth_texture.png>) | Earth texture used for the 3D plots |

**Result:** Strategy 3 (bi-elliptic) selected → **Δv = 3.2254 km/s**, ΔT = 2.062 × 10⁶ s (see [Table A.1–A.3](docs/report/Report_Space_LAB_PoliMi.pdf)).

🎬 Animation: [`docs/video/scenario1_animation.mp4`](docs/video/scenario1_animation.mp4)

### Scenario 2 — Interplanetary Direct Transfer to Nyx
`src/Project/Scenery 2/`

Direct heliocentric transfer between Earth's orbit and Nyx's orbit via a two-impulse elliptical transfer. The transfer geometry is fully described by three free parameters (departure true anomaly, arrival true anomaly, argument of pericenter of the transfer orbit); the optimal combination is found first with a 100×100×100 grid search, then refined with `fmincon` (SQP).

| File | Description |
|---|---|
| [`mission_transfer_nyx.m`](<src/Project/Scenery 2/mission_transfer_nyx.m>) | Main script: defines Earth/Nyx heliocentric orbits, sets up and solves the transfer-orbit optimization |
| [`grid_search_s2.m`](<src/Project/Scenery 2/grid_search_s2.m>) | 3-parameter grid search (θ₁, θ₂, ω_T) used as the initial guess for the `fmincon` refinement |
| [`scenery2_deltaV.m`](<src/Project/Scenery 2/scenery2_deltaV.m>) | Total Δv cost function (departure + arrival impulse) for a given transfer geometry |
| [`scenery2_plot.m`](<src/Project/Scenery 2/scenery2_plot.m>) | 3D plot of the Earth–Nyx transfer geometry (Fig. 2.2, 2.3) |
| [`sun_texture.jpg`](<src/Project/Scenery 2/sun_texture.jpg>) | Sun texture used for the 3D plots |

**Result:** optimal transfer orbit a_T ≈ 1.883 × 10⁸ km, e_T ≈ 0.2035 → **Δv = 7.0159 km/s** (departure 2.999 + arrival 4.017 km/s), **TOF = 61.33 days**.

🎬 Animation: [`docs/video/scenario2_animation.mp4`](docs/video/scenario2_animation.mp4)

### Scenario 3 — Escape and Arrival Phases
`src/Project/Scenery 3/`

Links Scenarios 1 and 2 via the **Patched Conics method**: a hyperbolic escape trajectory out of Earth's Sphere of Influence (analyzed both coplanar and non-coplanar) and a hyperbolic capture trajectory into Nyx's SOI.

| File | Description |
|---|---|
| [`scenery3_coplanar.m`](<src/Project/Scenery 3/scenery3_coplanar.m>) | Coplanar escape hyperbola (departure at parking-orbit pericenter) and Nyx capture hyperbola |
| [`scenery3_3D.m`](<src/Project/Scenery 3/scenery3_3D.m>) | Non-coplanar escape hyperbola: solves for the injection point/orientation matching the required v∞ direction, plus grid search over the optimal injection true anomaly |
| [`istogram_comparison_s3.m`](<src/Project/Scenery 3/istogram_comparison_s3.m>) | Bar-chart comparing pericenter-injection vs. optimal-injection Δv (Fig. 4.2) |
| [`scenery3_plot_coplanar.m`](<src/Project/Scenery 3/scenery3_plot_coplanar.m>) | 3D plot of the coplanar escape/capture hyperbolas (Fig. 3.2) |
| [`scenery3_plot_noncoplanar.m`](<src/Project/Scenery 3/scenery3_plot_noncoplanar.m>) | 3D plot of the non-coplanar escape hyperbola, optimal vs. pericenter injection (Fig. 3.1) |
| [`earth_texture.png`](<src/Project/Scenery 3/earth_texture.png>) | Earth texture used for the 3D plots |
| [`asteroid_texture.jpg`](<src/Project/Scenery 3/asteroid_texture.jpg>) | Nyx texture used for the 3D plots |

**Result:** Earth escape (optimal injection, θ = 165°) **Δv = 3.8304 km/s** vs. 15.09 km/s at the pericenter (**~74% saving**); Nyx capture **Δv = 4.0171 km/s** at r_p = 1.5 km.

🎬 Animation: [`docs/video/scenario3_optimal_injection.mp4`](docs/video/scenario3_optimal_injection.mp4)

---

## Documentation

| File | Content |
|---|---|
| [`docs/report/Report_Space_LAB_PoliMi.pdf`](docs/report/Report_Space_LAB_PoliMi.pdf) | Full final report: methodology, equations, results, and Appendix A with the complete orbital-parameter tables for every maneuver |
| [`docs/data/Scenery 1 data.pdf`](<docs/data/Scenery 1 data.pdf>), [`Scenery 2-3 data.pdf`](<docs/data/Scenery%202-3%20data.pdf>) | Group-specific numerical input data (orbital parameters, target body) assigned for the project |

## Methodology Summary

1. **Patched Conics Approximation** — the full N-body trajectory (Sun + Earth + Nyx) is split into three regions, each dominated by a single attractor (Earth's SOI, the Sun's SOI, Nyx's SOI), and solved as a sequence of independent two-body problems.
2. **Impulsive maneuver model** — every burn is modeled as an instantaneous velocity increment (Δv); no finite-burn or low-thrust effects are considered.
3. **Optimization** — each scenario's free design parameters (intermediate orbit shape, transfer geometry, injection point) are explored with a discrete grid search first, then, where applicable, refined with MATLAB's `fmincon` (Sequential Quadratic Programming) for higher precision.

### Known Limitations (see report, §4 for full discussion)
- All maneuvers are impulsive; no finite burn times.
- No orbital perturbations or mid-course correction burns are modeled.
- 3908 Nyx is modeled as a uniform-density point mass/sphere (no gravitational harmonics).
- The geocentric environment is idealized (no debris, no third-body perturbations).

---

## Getting Started

**Requirements:** MATLAB (developed with the base MATLAB environment; `fmincon` requires the **Optimization Toolbox**).

1. Clone the repository:
   ```bash
   git clone https://github.com/Bernaaaa/Mission-analysis-project.git
   cd Mission-analysis-project
   ```
2. Add the toolbox to your MATLAB path (or `cd` into it), since every project script depends on it:
   ```matlab
   addpath(genpath('src/lab'));
   ```
3. Run a scenario's main script, e.g.:
   ```matlab
   cd 'src/Project/Scenery 1'
   scenery1_standard      % Strategies 1 & 2
   scenery1_deltaV        % Strategy 3 (bi-elliptic) + grid-search optimization
   ```
   Repeat analogously for `Scenery 2` (`mission_transfer_nyx.m`) and `Scenery 3` (`scenery3_coplanar.m`, `scenery3_3D.m`).
4. Review the generated plots alongside [`docs/report/Report_Space_LAB_PoliMi.pdf`](docs/report/Report_Space_LAB_PoliMi.pdf) and the animations in [`docs/video/`](docs/video) for the full context and discussion of results.

---

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request with suggestions, bug fixes, or extensions (e.g. finite-burn modeling, N-body validation, alternative optimization strategies as suggested in the report's Conclusions).

## License

This project is licensed under the MIT License — see the [`license`](license) file for details.

## Contact

For questions about the project, please open an issue on the repository.
