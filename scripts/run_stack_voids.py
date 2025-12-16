#!/usr/bin/env python3
"""
run_stack_voids.py

Convenience script for stacking CMB temperatures around DESI voids
and reporting monopole/dipole diagnostics.

This script is intentionally lightweight and delegates all heavy
lifting to the gop_boundary_desi package.
"""

from pathlib import Path
import numpy as np
import healpy as hp

from gop_boundary_desi.io import read_void_catalog
from gop_boundary_desi.stacking import stack_cmb_temperature_profile
from gop_boundary_desi.dipole import fit_monopole_dipole


def main():
    # ---- USER PATHS (edit as needed) ----
    void_catalog = Path("data/desi_voids.fits")
    cmb_map_path = Path("data/planck_smica.fits")
    mask_path = Path("data/mask.fits")  # optional

    theta_max_deg = 5.0
    nbins = 15
    max_voids = None  # set to int for quick tests

    # ---- LOAD DATA ----
    print("Loading void catalog...")
    voids = read_void_catalog(void_catalog)

    print("Loading CMB map...")
    cmb_map = hp.read_map(cmb_map_path, field=0, verbose=False)
    nside = hp.get_nside(cmb_map)

    mask = None
    if mask_path.exists():
        print("Loading mask...")
        mask_map = hp.read_map(mask_path, field=0, verbose=False)
        mask = mask_map != 0

    # ---- STACKING ----
    theta_bins = np.linspace(0.0, theta_max_deg, nbins + 1)

    print("Stacking CMB temperatures around void centers...")
    result = stack_cmb_temperature_profile(
        voids=voids,
        cmb_map=cmb_map,
        nside=nside,
        theta_bins_deg=theta_bins,
        mask=mask,
        max_voids=max_voids,
    )

    print(f"\nUsed {result.n_used} voids")
    print("theta_deg   mean_uK   sem_uK")
    for th, mu, se in zip(
        result.theta_bins_deg,
        result.prof_mean_uK,
        result.prof_sem_uK,
    ):
        print(f"{th:8.3f} {mu:9.3f} {se:8.3f}")

    # ---- DIPOLE DIAGNOSTIC ----
    print("\nFitting monopole + dipole on masked CMB map...")
    dfit = fit_monopole_dipole(cmb_map, mask=mask)
    print(f"Dipole amplitude : {dfit.amp:.4e}")
    print(f"Dipole axis (x,y,z): {dfit.axis_vec}")
    print(f"chi2 / dof       : {dfit.chi2 / dfit.dof:.4e}")


if __name__ == "__main__":
    main()
