# GoP-Boundary-DESI-EUCLID

# Repo Layout
```
gop-boundary-desi/
  README.md
  requirements.txt
  gop_boundary_desi/
    io.py          # Input readers (void catalogs, maps)
    boundary.py    # Boundary / saturation templates
    stacking.py    # Void stacking routines
    dipole.py      # Monopole/dipole fitting and regression
    cli.py         # Command-line interface
  scripts/
    run_stack_voids.py
  notebooks/
    01_quickstart.ipynb

```

Installation

```
git clone https://github.com/your-username/gop-boundary-desi.git
cd gop-boundary-desi
pip install -r requirements.txt
```


#  Gravity of Probability Boundary Diagnostics for DESI DR2

This repository provides analysis tools to test whether decoherence-induced boundary saturation effects—as developed in the accompanying Gravity of Probability (GoP) framework—can explain large-scale dipole mismatches and void-stacking anomalies without introducing new microphysical degrees of freedom.

The code is designed as a diagnostic and validation layer, parallel to the theoretical paper and thesis, and is suitable for application to DESI DR2 preview products, Planck CMB maps, and future Euclid data.

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


This repository is designed to sit parallel to the theoretical work: 
(PLACEHOLDER - NOT YET PUBLISHED, POST-DESI VACs Validation of primary work
https://figshare.com/articles/thesis/The_Gravity_of_Probability_i_Replicating_Dark_Matter_Effects_Through_Quantum_Decoherence_Curvature_i_/29815934 

https://figshare.com/articles/preprint/DESI_DR2_VACs_Predictions/30593876?file=59479682 )

The paper establishes the mechanism: decoherence → curvature → boundary → fossilization

This code tests whether that mechanism leaves detectable imprints in real data

No claims in this repository exceed those made in the paper. All diagnostics are downstream of, and consistent with, the published framework.











