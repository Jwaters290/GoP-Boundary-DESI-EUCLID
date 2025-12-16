from __future__ import annotations

import numpy as np
import healpy as hp
from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass(frozen=True)
class DipoleFit:
    amp: float
    axis_vec: np.ndarray  # (3,)
    coeffs: np.ndarray    # (4,) [monopole, dx, dy, dz]
    chi2: float
    dof: int

def fit_monopole_dipole(map_in: np.ndarray, mask: Optional[np.ndarray] = None) -> DipoleFit:
    """
    Fit T(n) ≈ a0 + a·n (monopole + dipole) on unmasked pixels.
    """
    if mask is None:
        mask = np.ones_like(map_in, dtype=bool)
    m = np.asarray(map_in, dtype=float)
    good = mask & np.isfinite(m)

    nside = hp.get_nside(m)
    ipix = np.where(good)[0]
    vecs = np.array(hp.pix2vec(nside, ipix)).T  # (N, 3)

    X = np.column_stack([np.ones(vecs.shape[0]), vecs])  # (N,4)
    y = m[ipix]

    # Least squares
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ coeffs
    resid = y - yhat

    dof = max(y.size - 4, 1)
    chi2 = float(np.sum(resid**2))
    a = coeffs[1:]
    amp = float(np.linalg.norm(a))
    axis = a / amp if amp > 0 else np.array([np.nan, np.nan, np.nan])
    return DipoleFit(amp=amp, axis_vec=axis, coeffs=coeffs, chi2=chi2, dof=dof)

def regress_template(y: np.ndarray, t: np.ndarray, mask: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """
    Fit y ≈ b * t + c on masked pixels; returns (b, R^2).
    """
    if mask is None:
        mask = np.ones_like(y, dtype=bool)
    good = mask & np.isfinite(y) & np.isfinite(t)
    yy = y[good]
    tt = t[good]

    X = np.column_stack([tt, np.ones(tt.size)])
    coeffs, *_ = np.linalg.lstsq(X, yy, rcond=None)
    b, c = coeffs

    yhat = X @ coeffs
    ss_res = np.sum((yy - yhat) ** 2)
    ss_tot = np.sum((yy - np.mean(yy)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(b), float(r2)
