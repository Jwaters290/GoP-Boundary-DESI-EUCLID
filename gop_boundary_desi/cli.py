from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import healpy as hp

from .io import read_void_catalog
from .stacking import stack_cmb_temperature_profile
from .dipole import fit_monopole_dipole
from .boundary import BoundaryParams, make_boundary_template_map

def main() -> None:
    ap = argparse.ArgumentParser("gop-boundary-desi")
    ap.add_argument("--voids", required=True, help="Path to void catalog (CSV/Parquet/FITS)")
    ap.add_argument("--cmb", required=True, help="Path to HEALPix CMB map (FITS)")
    ap.add_argument("--mask", default=None, help="Optional HEALPix mask (FITS). Nonzero=keep.")
    ap.add_argument("--nside", type=int, default=None, help="Override NSIDE (else inferred from map)")
    ap.add_argument("--theta-max", type=float, default=5.0, help="Max stacking radius in deg")
    ap.add_argument("--nbins", type=int, default=15, help="Number of radial bins")
    ap.add_argument("--max-voids", type=int, default=None, help="Limit number of voids for quick tests")
    ap.add_argument("--template", action="store_true", help="Also build a boundary template and fit it")
    ap.add_argument("--rb", type=float, default=1.0, help="Boundary proxy radius")
    ap.add_argument("--dr", type=float, default=0.1, help="Boundary proxy thickness")
    ap.add_argument("--A", type=float, default=0.2, help="Dipole modulation amplitude")
    ap.add_argument("--axis-lon", type=float, default=264.0, help="Dipole axis lon (deg, Galactic)")
    ap.add_argument("--axis-lat", type=float, default=48.0, help="Dipole axis lat (deg, Galactic)")
    args = ap.parse_args()

    voids = read_void_catalog(args.voids)

    cmb = hp.read_map(args.cmb, field=0, verbose=False)
    nside = args.nside or hp.get_nside(cmb)

    mask = None
    if args.mask:
        m = hp.read_map(args.mask, field=0, verbose=False)
        mask = m != 0

    theta_bins = np.linspace(0.0, args.theta_max, args.nbins + 1)
    res = stack_cmb_temperature_profile(
        voids=voids,
        cmb_map=cmb,
        nside=nside,
        theta_bins_deg=theta_bins,
        mask=mask,
        max_voids=args.max_voids,
    )

    print(f"Stack used N={res.n_used} voids")
    print("theta_deg  mean_uK  sem_uK")
    for th, mu, se in zip(res.theta_bins_deg, res.prof_mean_uK, res.prof_sem_uK):
        print(f"{th:8.3f} {mu:9.3f} {se:8.3f}")

    # Dipole on the masked CMB map (baseline diagnostic)
    dfit = fit_monopole_dipole(cmb, mask=mask)
    print(f"\nCMB monopole+dipole fit: amp={dfit.amp:.4e} axis={dfit.axis_vec} chi2/dof={dfit.chi2/dfit.dof:.4e}")

    if args.template:
        # Convert (lon, lat) to unit vector in HEALPix convention (theta, phi)
        lon = np.deg2rad(args.axis_lon)
        lat = np.deg2rad(args.axis_lat)
        theta = np.pi / 2 - lat
        phi = lon
        n_hat = np.array(hp.ang2vec(theta, phi))

        # Example r_of_pix: proxy "depth" based on |sin(b)| to mimic survey boundary weighting.
        # Replace with a survey-specific selection-depth model or z-weighted effective r(z).
        def r_of_pix(v: np.ndarray) -> float:
            # v is unit vector. "Galactic latitude proxy"
            b = np.arcsin(v[2])
            return 1.0 + 0.5 * np.abs(np.sin(b))

        params = BoundaryParams(rb=args.rb, dr=args.dr, A=args.A, n_hat=n_hat)
        tmpl = make_boundary_template_map(nside=nside, params=params, r_of_pix=r_of_pix, normalize=True)

        # Save template for inspection
        out = Path("boundary_template.fits")
        hp.write_map(str(out), tmpl, overwrite=True)
        print(f"\nWrote boundary template: {out.resolve()}")

if __name__ == "__main__":
    main()
