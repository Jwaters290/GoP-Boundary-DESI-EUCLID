from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import healpy as hp

@dataclass(frozen=True)
class BoundaryParams:
    """
    Geometric/phenomenological boundary template parameters.
    These are NOT GoP microphysical parameters; they describe an observational
    boundary shell/gradient on the sky.

    rb: characteristic boundary radius (in comoving units or proxy units)
    dr: thickness/softening scale
    A: dipole modulation amplitude (dimensionless)
    n_hat: dipole axis as unit vector in Cartesian coordinates
    """
    rb: float
    dr: float
    A: float
    n_hat: np.ndarray  # shape (3,)

def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))

def heavy_boundary_weight(r: np.ndarray, rb: float, dr: float) -> np.ndarray:
    """
    Saturation shell/gradient weight as a function of radius r.
    For a 'heavy boundary', dr is small and the transition is sharp.
    """
    if dr <= 0:
        raise ValueError("dr must be > 0")
    return sigmoid((r - rb) / dr)

def dipole_modulation(n_vec: np.ndarray, n_hat: np.ndarray, A: float) -> np.ndarray:
    """
    1 + A cos(theta) where cos(theta)=n·n_hat.
    """
    n_hat = np.asarray(n_hat, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    cosang = np.clip(np.dot(n_vec, n_hat), -1.0, 1.0)
    return 1.0 + A * cosang

def make_boundary_template_map(
    nside: int,
    params: BoundaryParams,
    r_of_pix: callable,
    normalize: bool = True,
) -> np.ndarray:
    """
    Build a HEALPix map template T(n) = W(r(n)) * (1 + A n·n_hat).

    r_of_pix: user-supplied function mapping pixel direction vector -> effective r.
              In practice, you can use a redshift-weighted r(z) for the void sample
              or an analysis proxy (e.g., survey depth / selection boundary).
    """
    npix = hp.nside2npix(nside)
    vecs = np.array(hp.pix2vec(nside, np.arange(npix))).T  # (npix, 3)

    r = np.array([r_of_pix(v) for v in vecs], dtype=float)
    W = heavy_boundary_weight(r, params.rb, params.dr)
    M = np.array([dipole_modulation(v, params.n_hat, params.A) for v in vecs], dtype=float)
    temp = W * M

    if normalize:
        temp = temp - np.mean(temp)
        s = np.std(temp)
        if s > 0:
            temp = temp / s
    return temp
