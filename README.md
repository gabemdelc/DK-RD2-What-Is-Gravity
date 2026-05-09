# DK-RD2 — What Is Gravity?

**From Emergent Mass to Thermodynamic–Geometric Dynamics in the DK-RD2 Framework**

Author: **Gabriel Martín del Campo Flores**  
Signature: **GabE=mc²**  
Zenodo DOI: https://doi.org/10.5281/zenodo.20078175  
GitHub: https://github.com/gabemdelc/DK-RD2-What-Is-Gravity

---

## Overview

This repository contains the reproducible Python code used for the paper:

**“What Is Gravity? From Emergent Mass to Thermodynamic–Geometric Dynamics in the DK-RD2 Framework”**

The DK-RD2 framework explores gravity as an emergent thermodynamic–geometric interaction rather than a fixed fundamental coupling. In this model, the effective gravitational coupling depends on the thermodynamic and relativistic state of the system:

```text
Gab(T, v) ≈ G0 · [1 + αDK · (v²/c²) · (T0/T)]
````

where:

* `G0` is Newton’s gravitational constant.
* `T0` is the present-day CMB temperature.
* `T` is the effective thermodynamic state.
* `v` is the characteristic relativistic velocity scale.
* `αDK` is the effective geometric projection factor calibrated from observational data.

The model investigates whether cosmic acceleration, emergent de Sitter behavior, compact gravitational systems, and effective dark-energy-like signatures can arise from thermodynamic–geometric evolution without introducing a fundamental cosmological constant.

---

## Scientific Goals

This codebase implements the numerical and observational tests associated with DK-RD2:

* DESI BAO + Cosmic Chronometer calibration of `αDK`
* Statistical profiling of the DK projection factor
* Thermodynamic–relativistic gravitational coupling curves
* Effective equation-of-state reconstruction `w_eff(z)`
* Thermodynamic activity versus projected gravitational coupling
* Compact SIDM-halo effective-gravity reinterpretation
* Reproducible figure and CSV evidence generation

---

## Main Physical Idea

DK-RD2 proposes that gravity emerges from the organization of energy across thermodynamic, relativistic, and geometric configurations.

In this interpretation:

```text
Energy Distribution → Localization → Interaction → Saturation → Evaporation → Degravification
```

Gravity is treated as a projected macroscopic interaction arising from the thermodynamic organization of energy, rather than as a primitive fixed-force field.

---

## Repository Structure

Recommended structure:

```text
DK-RD2-What-Is-Gravity/
│
├── DK_RD2_Core.py
├── DK-RD2_ThermoGravity_Framework.py
├── README.md
├── requirements.txt
│
├── data/
│   └── DESI/
│       └── bao_data/
│           └── desi_bao_dr2/
│               ├── desi_gaussian_bao_ALL_GCcomb_mean.txt
│               └── desi_gaussian_bao_ALL_GCcomb_cov.txt
│
├── evidence/
│   ├── figures/
│   ├── tables/
│   └── stats/
│
└── paper/
    └── What_Is_Gravity_DK-RD2_Framework.pdf
```

---

## Main Files

### `DK_RD2_Core.py`

Core computational engine.

Includes:

* Physical constants
* Thermodynamic gravitational coupling `Gab(T, v)`
* ΛCDM reference functions
* DK-RD2 expansion functions
* Distance calculations
* BAO utilities
* Lensing utilities
* Effective equation-of-state reconstruction
* DESI BAO + CC calibration helpers

This file is designed to remain mostly model/core logic only.

---

### `DK-RD2_ThermoGravity_Framework.py`

Main driver script for the paper.

Includes:

* Dataset loading
* Figure generation
* Statistical calibration
* Evidence table export
* CSV/statistical output generation
* Paper-specific plotting configuration

This is the main executable script used to reproduce the figures and evidence tables.

---

## Data Requirements

This project uses processed DESI BAO DR2 likelihood data.

Expected files:

```python
path = "data/DESI/bao_data/desi_bao_dr2/"

desi_bao_mean_path = path + "desi_gaussian_bao_ALL_GCcomb_mean.txt"
desi_bao_cov_path  = path + "desi_gaussian_bao_ALL_GCcomb_cov.txt"
```

These files correspond to the public DESI BAO likelihood compilation distributed through the CobayaSampler BAO data repository.

No CLASS or Boltzmann solver is required for the current version of this paper.

---

## Installation

Create a Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install dependencies:

```bash
pip install numpy pandas scipy matplotlib
```

Optional:

```bash
pip install astropy
```

---

## Running the Code

From the repository root:

```bash
python DK-RD2_ThermoGravity_Framework.py
```

The script generates:

* figures
* CSV tables
* statistical summaries
* reproducibility metadata

Outputs are written under:

```text
evidence/
```

---

## Figures Generated

### Figure 01 — DESI Expansion and BAO Calibration of αDK

Joint DESI BAO + Cosmic Chronometer calibration of the DK-RD2 projection factor `αDK` and the sound horizon scale `rd`.

---

### Figure 02 — Statistical Constraint on the DK Projection Factor

Profiled `χ²(αDK)` likelihood curve with `rd` independently profiled for each value of `αDK`.

This figure demonstrates that `αDK` is observationally constrained rather than arbitrary.

---

### Figure 03 — Thermodynamic–Geometric Structure of the DK-RD2 Coupling

Comparison between:

```text
Gab_max / G0
Gab_eff / G0
```

showing bounded saturation and recovery of the Newtonian limit.

---

### Figure 04 — Diagnostic Reconstruction of w_eff(z)

Effective FRW-style reconstruction of the equation of state that would be inferred if DK-RD2 expansion were interpreted as dark energy.

This does not imply a fundamental dark-energy fluid.

---

### Figure 05 — Thermodynamic Activity vs Projected Gravitational Coupling

Comparison between energetic confinement, thermal activity, and projected gravitational coupling.

Includes consistency checks for the projection identity.

---

### Figure 06 — Compact SIDM Halo Effective-Gravity Reinterpretation

Comparison between published compact SIDM halo profiles and an effective-gravity reinterpretation within DK-RD2.

The original SIDM interpretation is preserved. DK-RD2 is introduced only as an alternative effective gravitational representation.

---

## Reproducibility

The framework exports numerical evidence tables for each figure.

Typical outputs include:

```text
evidence/figure_XX.png
evidence/table_XX.csv
evidence/table_XX_stats.csv
```

These files are intended to make the model directly testable, reproducible, and independently verifiable.

---

## Scientific Status

DK-RD2 should be understood as an effective thermodynamic–geometric framework under development.

It is not presented as a completed fundamental theory.

Open areas include:

* first-principles derivation of `αDK`
* covariant relativistic formulation
* microscopic interpretation of layered spacetime
* structure formation predictions
* compact-object dynamics
* quantum-gravitational connection

---

## Citation

If you use this code or discuss the framework, please cite:

```text
Martín del Campo Flores, G. (2026).
What Is Gravity? From Emergent Mass to Thermodynamic–Geometric Dynamics
in the DK-RD2 Framework.
Zenodo.
https://doi.org/10.5281/zenodo.20078175
```

---

## License

MIT License.

---

## Disclaimer

This repository contains an independent theoretical and computational research framework.

The results are intended for scientific discussion, reproducibility, and falsifiable testing. Inclusion in public repositories or Zenodo communities does not imply peer review or endorsement of the claims.

```

Está alineado con el paper: la formulación de `Gab(T,v)` y `αDK`, la calibración DESI BAO + CC, las seis figuras y la aclaración de que el modelo sigue incompleto aparecen en el documento. :contentReference[oaicite:0]{index=0}
```
