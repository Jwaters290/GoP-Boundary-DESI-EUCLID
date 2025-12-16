#!/usr/bin/env python3
"""
GoP diagnostics: P(k) modulation and DeltaSigma(R) slope diagnostics.

This is a practice diagnostic script: it does not tune parameters or claim detection; it outputs “what GoP would do” under a clearly defined mapping from (k or R) → an effective decoherence energy argument.

This script is designed to be:
- runnable with no external datasets (uses analytic baselines)
- able to consume baseline data if provided (CSV)
- parameter-rigid: core GoP parameters are inputs; no fitting is performed

Outputs:
- CSV with k, P0(k), modulation M(k), P_gop(k)
- CSV with R, DS0(R), modulation M(R), DS_gop(R)
- Printed diagnostic summaries (amplitude metrics, slope changes)

Usage examples:
  # P(k) diagnostic (analytic baseline)
  python scripts/gop_diagnostics.py pk \
    --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
    --kmin 1e-3 --kmax 1 --nk 200 --out pk_diag.csv

  # P(k) diagnostic with baseline CSV (columns: k, P0)
  python scripts/gop_diagnostics.py pk \
    --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
    --pk-baseline baseline_pk.csv --out pk_diag.csv

  # DeltaSigma diagnostic (analytic baseline)
  python scripts/gop_diagnostics.py ds \
    --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
    --Rmin 0.05 --Rmax 50 --nR 200 --out ds_diag.csv

  # DeltaSigma with baseline CSV (columns: R, DS0)
  python scripts/gop_diagnostics.py ds \
    --kappaA 1.5e-15 --E0 1e12 --f_ent 0.20 --A_CP 0.0245 \
    --ds-baseline baseline_ds.csv --out ds_diag.csv
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


# -----------------------------
# Core GoP kernel + diagnostics
# -----------------------------

@dataclass(frozen=True)
class GoPParams:
    kappaA: float   # s^-1 erg^-1 (or your chosen consistent units)
    E0: float       # erg
    f_ent: float    # dimensionless
    A_CP: float     # dimensionless (small)


def gamma_kernel(E: np.ndarray, p: GoPParams) -> np.ndarray:
    """
    Gamma(E) = kappaA * E * exp(1 - E/E0)
    """
    E = np.asarray(E, dtype=float)
    # Guard against negative/zero energies in any proxy mapping
    E = np.maximum(E, 1e-300)
    return p.kappaA * E * np.exp(1.0 - (E / p.E0))


def normalized_gamma(E: np.ndarray, p: GoPParams) -> np.ndarray:
    """
    Normalizes Gamma(E) by its maximum value over E, which occurs at E = E0:
      Gamma_max = kappaA * E0 * exp(1 - 1) = kappaA * E0
    So Gamma_norm(E) = Gamma(E) / (kappaA * E0)
    """
    g = gamma_kernel(E, p)
    gmax = p.kappaA * p.E0
    if gmax <= 0:
        raise ValueError("kappaA*E0 must be positive for normalization.")
    return g / gmax


def cp_modulation_factor(mu: np.ndarray, A_CP: float) -> np.ndarray:
    """
    A minimal CP-asymmetry modulation that can encode a parity-odd preference
    without introducing new microphysics.

    mu is cos(theta) relative to a chosen axis. For isotropic diagnostics, set mu=0.
    """
    mu = np.asarray(mu, dtype=float)
    # Small, bounded, sign-sensitive factor
    return 1.0 + A_CP * mu


def pk_modulation(
    k: np.ndarray,
    p: GoPParams,
    *,
    k0: float = 0.2,
    alpha_k: float = 1.0,
    mu: float = 0.0,
    strength: float = 1.0
) -> np.ndarray:
    """
    Returns multiplicative modulation M(k) applied to baseline P0(k).

    We define an energy proxy:
      E(k) = E0 * (k/k0)^alpha_k

    Then:
      M(k) = 1 + strength * f_ent * Gamma_norm(E(k)) * CP(mu)

    Notes:
    - This is a diagnostic template; k0 and alpha_k should be documented choices.
    - mu=0 yields isotropic CP-neutral diagnostic; nonzero mu can probe dipole-like modulation.
    """
    k = np.asarray(k, dtype=float)
    k = np.maximum(k, 1e-300)
    E = p.E0 * (k / k0) ** alpha_k
    gN = normalized_gamma(E, p)
    cp = cp_modulation_factor(np.array(mu), p.A_CP)
    return 1.0 + strength * p.f_ent * gN * cp


def ds_modulation(
    R: np.ndarray,
    p: GoPParams,
    *,
    R0: float = 1.0,
    alpha_R: float = 1.0,
    mu: float = 0.0,
    strength: float = 1.0
) -> np.ndarray:
    """
    Returns multiplicative modulation M(R) applied to baseline DeltaSigma0(R).

    Energy proxy:
      E(R) = E0 * (R0/R)^alpha_R

    Then:
      M(R) = 1 + strength * f_ent * Gamma_norm(E(R)) * CP(mu)

    Diagnostic intent: check slope changes or flattening/steepening across R.
    """
    R = np.asarray(R, dtype=float)
    R = np.maximum(R, 1e-300)
    E = p.E0 * (R0 / R) ** alpha_R
    gN = normalized_gamma(E, p)
    cp = cp_modulation_factor(np.array(mu), p.A_CP)
    return 1.0 + strength * p.f_ent * gN * cp


# -----------------------------
# Baseline loaders (CSV)
# -----------------------------

def load_two_column_csv(path: str, xcol: str, ycol: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Loads a CSV with headers; returns x, y.
    Expected headers: xcol, ycol.
    """
    import csv
    xs, ys = [], []
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f)
        if xcol not in r.fieldnames or ycol not in r.fieldnames:
            raise ValueError(f"CSV must contain columns '{xcol}' and '{ycol}'. Found: {r.fieldnames}")
        for row in r:
            xs.append(float(row[xcol]))
            ys.append(float(row[ycol]))
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    # sort by x
    s = np.argsort(x)
    return x[s], y[s]


def write_csv(path: str, header: list[str], cols: list[np.ndarray]) -> None:
    import csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i in range(len(cols[0])):
            w.writerow([float(c[i]) for c in cols])


# -----------------------------
# Diagnostics
# -----------------------------

def diagnostic_pk(k: np.ndarray, P0: np.ndarray, M: np.ndarray) -> dict:
    """
    Returns robust scalar diagnostics for P(k) modulation.
    """
    Pg = P0 * M
    frac = (Pg - P0) / np.maximum(P0, 1e-300)  # fractional change
    # summarize over range
    out = {
        "frac_mean": float(np.mean(frac)),
        "frac_rms": float(np.sqrt(np.mean(frac**2))),
        "frac_min": float(np.min(frac)),
        "frac_max": float(np.max(frac)),
        "M_min": float(np.min(M)),
        "M_max": float(np.max(M)),
    }
    return out


def log_slope(x: np.ndarray, y: np.ndarray, x1: float, x2: float) -> float:
    """
    Compute slope d ln y / d ln x over [x1, x2] using linear regression in log-log.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = (x >= x1) & (x <= x2) & (y > 0) & np.isfinite(y) & np.isfinite(x)
    if np.sum(m) < 5:
        raise ValueError("Not enough valid points in the requested range to compute slope.")
    X = np.log(x[m])
    Y = np.log(y[m])
    A = np.vstack([X, np.ones_like(X)]).T
    slope, intercept = np.linalg.lstsq(A, Y, rcond=None)[0]
    return float(slope)


def diagnostic_ds(R: np.ndarray, DS0: np.ndarray, M: np.ndarray, R_lo: float, R_hi: float) -> dict:
    """
    Returns slope diagnostics for DeltaSigma(R).
    """
    DSg = DS0 * M
    s0 = log_slope(R, DS0, R_lo, R_hi)
    sg = log_slope(R, DSg, R_lo, R_hi)
    out = {
        "slope_baseline": s0,
        "slope_gop": sg,
        "delta_slope": sg - s0,
        "M_min": float(np.min(M)),
        "M_max": float(np.max(M)),
    }
    return out


# -----------------------------
# CLI
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="GoP diagnostics for P(k) modulation and DeltaSigma(R) slope.")
    sub = p.add_subparsers(dest="mode", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--kappaA", type=float, required=True, help="GoP kappaA (units consistent with E0)")
    common.add_argument("--E0", type=float, required=True, help="GoP E0 (erg)")
    common.add_argument("--f_ent", type=float, required=True, help="GoP entanglement fraction")
    common.add_argument("--A_CP", type=float, required=True, help="GoP CP asymmetry parameter")
    common.add_argument("--mu", type=float, default=0.0, help="cos(theta) relative to chosen axis (default 0)")
    common.add_argument("--strength", type=float, default=1.0, help="diagnostic scaling for template visibility (default 1)")

    # P(k)
    pk = sub.add_parser("pk", parents=[common], help="P(k) modulation diagnostic")
    pk.add_argument("--pk-baseline", type=str, default=None, help="CSV with columns: k,P0")
    pk.add_argument("--kmin", type=float, default=1e-3)
    pk.add_argument("--kmax", type=float, default=1.0)
    pk.add_argument("--nk", type=int, default=200)
    pk.add_argument("--k0", type=float, default=0.2, help="pivot k0 for E(k) proxy")
    pk.add_argument("--alpha-k", type=float, default=1.0, help="power alpha_k for E(k)=E0*(k/k0)^alpha_k")
    pk.add_argument("--ns", type=float, default=0.965, help="analytic baseline tilt if no baseline file")
    pk.add_argument("--out", type=str, default="pk_diag.csv")

    # DeltaSigma
    ds = sub.add_parser("ds", parents=[common], help="DeltaSigma(R) slope diagnostic")
    ds.add_argument("--ds-baseline", type=str, default=None, help="CSV with columns: R,DS0")
    ds.add_argument("--Rmin", type=float, default=0.05)
    ds.add_argument("--Rmax", type=float, default=50.0)
    ds.add_argument("--nR", type=int, default=200)
    ds.add_argument("--R0", type=float, default=1.0, help="pivot R0 for E(R) proxy")
    ds.add_argument("--alpha-R", type=float, default=1.0, help="power alpha_R for E(R)=E0*(R0/R)^alpha_R")
    ds.add_argument("--eta", type=float, default=0.8, help="analytic baseline DS0 ~ R^-eta if no baseline file")
    ds.add_argument("--slope-range", type=str, default="0.2,5.0", help="R_lo,R_hi range for slope diagnostic")
    ds.add_argument("--out", type=str, default="ds_diag.csv")

    return p


def main() -> None:
    args = build_parser().parse_args()
    gop = GoPParams(kappaA=args.kappaA, E0=args.E0, f_ent=args.f_ent, A_CP=args.A_CP)

    if args.mode == "pk":
        if args.pk_baseline:
            k, P0 = load_two_column_csv(args.pk_baseline, "k", "P0")
        else:
            k = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)
            # Analytic baseline (placeholder): P0(k) ~ k^(ns) * exp(-(k/kd)^2)
            # The exponential cutoff prevents blow-ups at high k in this simple diagnostic.
            kd = 2.0
            P0 = (k ** args.ns) * np.exp(-(k / kd) ** 2)

        M = pk_modulation(
            k, gop,
            k0=args.k0,
            alpha_k=args.alpha_k,
            mu=args.mu,
            strength=args.strength
        )
        Pg = P0 * M

        d = diagnostic_pk(k, P0, M)

        print("P(k) modulation diagnostic")
        print("-------------------------")
        print(f"kappaA       : {gop.kappaA:.6e}")
        print(f"E0           : {gop.E0:.6e}")
        print(f"f_ent        : {gop.f_ent:.6g}")
        print(f"A_CP         : {gop.A_CP:.6g}")
        print(f"mu           : {args.mu:.6g}")
        print(f"k0, alpha_k  : {args.k0:.6g}, {args.alpha_k:.6g}")
        print("")
        print(f"M_min, M_max : {d['M_min']:.6g}, {d['M_max']:.6g}")
        print(f"frac_mean    : {d['frac_mean']:.6g}")
        print(f"frac_rms     : {d['frac_rms']:.6g}")
        print(f"frac_min/max : {d['frac_min']:.6g}, {d['frac_max']:.6g}")
        print("")
        print(f"Writing CSV: {args.out}")

        write_csv(
            args.out,
            header=["k", "P0", "M", "P_gop"],
            cols=[k, P0, M, Pg]
        )

    elif args.mode == "ds":
        if args.ds_baseline:
            R, DS0 = load_two_column_csv(args.ds_baseline, "R", "DS0")
        else:
            R = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nR)
            # Analytic baseline (placeholder): DS0(R) ~ R^-eta
            DS0 = R ** (-args.eta)

        M = ds_modulation(
            R, gop,
            R0=args.R0,
            alpha_R=args.alpha_R,
            mu=args.mu,
            strength=args.strength
        )
        DSg = DS0 * M

        R_lo, R_hi = [float(x) for x in args.slope_range.split(",")]
        d = diagnostic_ds(R, DS0, M, R_lo=R_lo, R_hi=R_hi)

        print("DeltaSigma(R) slope diagnostic")
        print("------------------------------")
        print(f"kappaA       : {gop.kappaA:.6e}")
        print(f"E0           : {gop.E0:.6e}")
        print(f"f_ent        : {gop.f_ent:.6g}")
        print(f"A_CP         : {gop.A_CP:.6g}")
        print(f"mu           : {args.mu:.6g}")
        print(f"R0, alpha_R  : {args.R0:.6g}, {args.alpha_R:.6g}")
        print(f"slope range  : {R_lo:g} to {R_hi:g}")
        print("")
        print(f"M_min, M_max : {d['M_min']:.6g}, {d['M_max']:.6g}")
        print(f"slope_baseline : {d['slope_baseline']:.6g}")
        print(f"slope_gop      : {d['slope_gop']:.6g}")
        print(f"delta_slope    : {d['delta_slope']:.6g}")
        print("")
        print(f"Writing CSV: {args.out}")

        write_csv(
            args.out,
            header=["R", "DS0", "M", "DS_gop"],
            cols=[R, DS0, M, DSg]
        )

    else:
        raise RuntimeError("Unknown mode.")


if __name__ == "__main__":
    main()
