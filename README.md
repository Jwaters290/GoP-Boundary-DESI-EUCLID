# Parallel Repositories
All repositories published under the same account operate in parallel to this repository and are grounded in the same core theoretical claims of the Gravity of Probability (GoP) framework.

Each repository explores a distinct observational, numerical, or methodological aspect of the theory, but no repository introduces independent or conflicting assumptions. Results, diagnostics, and analyses are intended to be mutually consistent and cross-verifiable.

A full list of related repositories can be found at:

https://github.com/Jwaters290

This repository focuses specifically on boundary diagnostics and large-scale anisotropy tests using DESI and EUCLID data products, and should be interpreted as a computational companion to the primary Gravity of Probability papers and thesis.

## Citation

If you use this code or methodology, please cite:

Waters, J. (2025). *Gravity of Probability: Boundary Diagnostics Using DESI and EUCLID*.  
GitHub repository: https://github.com/Jwaters290/GoP-Boundary-DESI-EUCLID



# GoP(Gravity of Probability)-Boundary-DESI-EUCLID

# Repo Layout
```
GoP-Boundary-DESI-EUCLID/
├─ .gitignore
├─ LICENSE
├─ README.md                  # Detailed overview, motivation, usage, data requirements
├─ requirements.txt           # (Referenced; likely numpy, scipy, astropy, healpy, pandas, matplotlib, tqdm)
├─ gop_boundary_desi/         # Core package
│   ├─ __init__.py
│   ├─ io.py                  # Void catalog/CMB map readers (CSV/Parquet/FITS/HEALPix)
│   ├─ boundary.py            # Boundary/saturation template builders (shells/gradients)
│   ├─ stacking.py            # Void stacking for CMB ΔT profiles
│   ├─ dipole.py              # Monopole/dipole fits + template regression
│   └─ cli.py                 # Command-line interface orchestrator
├─ scripts/
│   └─ run_stack_voids.py     # Example workflow script
└─ notebooks/
    └─ 01_quickstart.ipynb    # Introductory walkthrough (install, data prep, basic run)

```

Installation

```
git clone [https://github.com/your-username/gop-boundary-desi.git](https://github.com/Jwaters290/GoP-Boundary-DESI-EUCLID)
cd gop-boundary-desi
pip install -r requirements.txt
```
# Implementation Notes and Future Extensions

The current implementation prioritizes clarity, transparency, and conservative assumptions over optimization or model complexity. Several extensions are intentionally left inactive or unimplemented to avoid introducing additional degrees of freedom prior to validated data releases.

Planned or optional extensions include:

 - Performance optimizations for very large void catalogs (e.g., vectorized HEALPix queries or cached ring selections).
 - Alternative boundary depth proxies (e.g., redshift-weighted comoving distance or survey depth maps), with the current latitude-based proxy retained as a neutral geometric placeholder.
 - Optional visualization outputs (e.g., saved profile plots or regression diagnostics), disabled by default to avoid presentation bias.
 - Synthetic validation tests using mock void catalogs and injected dipole signals to verify recovery behavior under controlled conditions.
 - These extensions are not required for the core diagnostic goals of the repository and will be introduced only when they improve robustness without compromising falsifiability.


# GoP-Boundary-DESI-EUCLID/scripts/gop_practice_diagnostics.py

Swap in real data (no code changes)

This script uses analytic placeholder baselines unless you provide a baseline CSV. To use real data later, just:

 - export a baseline CSV
 - pass it into the script
No GoP parameters or code changes are required.


A) Use a real baseline for P(k)

Baseline CSV format (headers required):

 - k, P0

```
python scripts/gop_practice_diagnostics.py pk \
  --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
  --pk-baseline data/baselines/pk_baseline.csv \
  --out outputs/pk_diag.csv
```

B) Use a real baseline for ΔΣ(R)

Baseline CSV format (headers required):

 - R, DS0

```
python scripts/gop_practice_diagnostics.py ds \
  --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
  --ds-baseline data/baselines/deltasigma_baseline.csv \
  --out outputs/ds_diag.csv
```



#  Gravity of Probability Boundary Diagnostics for DESI DR2

This repository provides analysis tools to test whether decoherence-induced boundary saturation effects—as developed in the accompanying Gravity of Probability (GoP) framework—can explain large-scale dipole mismatches and void-stacking anomalies without introducing new microphysical degrees of freedom.

The code is designed as a diagnostic and validation layer, parallel to the theoretical paper and thesis, and is suitable for application to DESI DR2 preview products, Planck CMB maps, and future Euclid data.

DESI Data Release 2 (DR2) cosmology results were formally published on March 19, 2025, including measurements of baryon acoustic oscillations (BAO) from both the Lyman-alpha forest and combined galaxy/quasar samples, and associated cosmological constraints. 
DESI Data

Supporting validation and methodological papers for these results are also available.

At the time of writing, specific value-added catalogs (VACs) targeting void stacking, lensing cross-products, or large-angle residual diagnostics have not yet been released.
Value-added Catalogs: https://data.desi.lbl.gov/doc/releases/edr/#value-added-catalogs

When such VACs become public, they can be incorporated into this repository’s workflows by providing appropriate baseline CSVs (e.g., k,P0 for power spectra or R,DS0 for lensing profiles) to the diagnostic scripts.

See also: DESI DR2 Publications — https://data.desi.lbl.gov/doc/papers/dr2/


# Validation Philosophy
This repository is designed as a diagnostic and falsifiable analysis pipeline, not a parameter-tuning or model-fitting framework.

All Gravity of Probability (GoP) parameters used here are globally fixed and derived outside this codebase. The boundary and saturation templates implemented in this repository are phenomenological diagnostics only, intended to test whether specific large-scale geometric imprints—if present—are consistent with decoherence-driven boundary effects predicted by the GoP framework.

No claims in this repository depend on positive detections.
Both null results and negative regressions are scientifically meaningful and serve to constrain or falsify boundary interpretations. Any apparent reduction in residuals or improved regression metrics must be reproducible across datasets and analysis choices to be considered physically significant.

This separation between theory, prediction, and diagnostic testing is deliberate and intended to minimize confirmation bias.

# Scientific Motivation

Recent experimental and theoretical work has established decoherence as a physical, time-directed process that converts quantum possibility into classical outcomes. Within the Gravity of Probability framework, this process sources a conserved probabilistic stress-energy contribution to spacetime curvature.

A key result of the framework is that decoherence acts as a boundary-setting process: once probabilistic curvature is generated, it remains fossilized in the spacetime metric even after decoherence saturates. Geometry therefore retains memory of where and how decoherence completed.

At cosmological scales, this implies that:
 - large-scale survey boundaries,
 - low-density environments (cosmic voids),
 - and anisotropic decoherence histories

can imprint persistent, direction-dependent curvature signatures that manifest as dipole residuals, void temperature anomalies, or lensing offsets—without invoking nonbaryonic dark matter or additional fundamental physics.

This repository implements analysis tools to test that hypothesis directly against observational data.

# Scope and Philosophy

This codebase is intentionally conservative in scope.

 - It does not introduce or tune any new microphysical parameters.
 - It does not modify the underlying GoP kernel or stress-energy tensor.
 - It treats boundary effects as geometric and observational diagnostics, not as new theory inputs.

The four core GoP quantities
(kappaA, E0, entanglement fraction, CP asymmetry)
are assumed fixed by first-principles physics and are not adjusted here.

All additional quantities used in this repository (e.g., boundary shells, gradients, dipole templates) are derived or phenomenological descriptors used only to test whether boundary saturation can account for observed anomalies.

# What This Code Does
 - The repository provides tools to:
 - Load DESI void catalogs (CSV, Parquet, or FITS)
 - Load HEALPix CMB temperature maps (e.g., Planck SMICA)
 - Stack CMB temperature profiles around void centers
 - Measure monopole and dipole residuals on masked skies
 - Construct “heavy boundary” or saturation-gradient templates
 - Regress observed residuals against boundary templates
 - Quantify whether boundary geometry explains dipole mismatches

This allows a direct test of the claim:

A strong decoherence boundary or saturation shell can bias large-angle observables and void stacks in a way consistent with data, without introducing new physics.

# Data Inputs
You will need:

1. A void catalog from DESI (DR2 preview or later)
 - Required columns: RA (deg), DEC (deg), Z
 - Optional: void radius (Mpc/h)

2. A CMB temperature map in HEALPix format
 - Planck SMICA is recommended for baseline tests

3. An optional sky mask
 - Galactic mask, point source mask, or DESI footprint

This repository does not redistribute DESI or Planck data. Users must obtain data from official sources.

#  Example Usage

Basic void stacking and dipole diagnostics:

```
python -m gop_boundary_desi.cli \
  --voids /path/to/desi_voids.fits \
  --cmb /path/to/planck_smica.fits \
  --mask /path/to/mask.fits \
  --theta-max 5 \
  --nbins 15
```
Including construction and fitting of a boundary template:
```
python -m gop_boundary_desi.cli \
  --voids /path/to/desi_voids.fits \
  --cmb /path/to/planck_smica.fits \
  --mask /path/to/mask.fits \
  --theta-max 5 \
  --nbins 15 \
  --template \
  --rb 1.0 \
  --dr 0.1 \
  --A 0.2
```
These parameters control only the shape of the diagnostic boundary template, not any fundamental GoP inputs.

# Interpretation of Results

Positive results (e.g., improved dipole alignment or reduced residuals after boundary regression) should be interpreted as:

 - Evidence that boundary saturation effects are observationally relevant
 - Support for the decoherence-boundary narrative developed in the paper
 - Motivation for more detailed survey-specific modeling

Negative results are equally valuable and would indicate that boundary geometry alone is insufficient, thereby sharpening falsifiability.


# Relationship to the Gravity of Probability Papers and Thesis


This repository is designed to sit parallel to the Gravity of Probability theoretical program rather than to replace or extend it.

The underlying framework—deriving probabilistic curvature from quantum decoherence and its boundary saturation—has been developed in a series of preprints and a publicly archived thesis, including:
 - The Gravity of Probability: Replicating Dark Matter Effects Through Quantum Decoherence Curvature (Figshare thesis archive)
https://figshare.com/articles/thesis/The_Gravity_of_Probability_i_Replicating_Dark_Matter_Effects_Through_Quantum_Decoherence_Curvature_i_/29815934
 - DESI DR2 VACs Predictions (pre-registered observational forecasts)
https://figshare.com/articles/preprint/DESI_DR2_VACs_Predictions/30593876



Those works establish the physical mechanism:

decoherence → probabilistic curvature → geometric boundary saturation → fossilized large-scale signatures

The purpose of this repository is narrower and complementary. It provides diagnostic tools to test whether that mechanism leaves detectable imprints in real survey data, such as void-stacked CMB temperature profiles and large-angle dipole residuals.

No claims in this repository exceed those made in the theoretical work. All diagnostics implemented here are downstream of, and fully consistent with, the published Gravity of Probability framework.



# Keywords

Gravity of Probability
Probabilistic curvature
Quantum decoherence gravity
Decoherence-induced curvature
Quantum information and gravity
Semiclassical gravity extension
Probabilistic stress-energy tensor
Interference-weighted quantum dynamics

Dark matter alternatives
No dark matter cosmology
Emergent gravity diagnostics
Non-particle dark matter
Information-based gravity
Decoherence-driven gravity

Cosmic voids
Void stacking analysis
DESI DR2
DESI void catalog
DESI VACs
Large-scale structure
Low-density cosmology
Cosmic boundary effects

CMB temperature anomalies
CMB dipole residuals
CMB large-angle anomalies
Planck SMICA
HEALPix analysis
Spherical harmonic analysis

Gravitational lensing diagnostics
Weak lensing anomalies
Boundary saturation effects
Geometric fossilization
Curvature boundary shell
Decoherence boundary theory

Pre-registered predictions
Falsifiable cosmology
Model-independent diagnostics
Survey boundary effects
Cosmological systematics testing

Python cosmology
Astrophysics data analysis
Healpy
NumPy SciPy cosmology
Open science cosmology
Reproducible research




