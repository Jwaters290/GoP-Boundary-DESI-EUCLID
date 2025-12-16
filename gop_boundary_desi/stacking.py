from __future__ import annotations

import numpy as np
import healpy as hp
from dataclasses import dataclass
from typing import Optional, Tuple

from .io import VoidCatalog

@dataclass(frozen=True)
class StackResult:
    theta_bins_deg: np.ndarray
    prof_mean_uK: np.ndarray
    prof_sem_uK: np.ndarray
    n_used: int

def _angdist_rad(v1: np.ndarray, v2: np.ndarray) -> float:
    return float(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

def stack_cmb_temperature_profile(
    voids: VoidCatalog,
    cmb_map: np.ndarray,
    nside: int,
    theta_bins_deg: np.ndarray,
    mask: Optional[np.ndarray] = None,
    uK_scale: float = 1e6,
    max_voids: Optional[int] = None,
) -> StackResult:
    """
    Stack CMB temperature around void centers.
    Uses ring-averaging by sampling pixels within theta bins around each void.

    NOTE: For speed/scaling, you’ll likely replace this with:
      - precomputed disc/ring pixel lists,
      - vectorized query_disc,
      - or a harmonic-space matched filter.
    """
    if mask is None:
        mask = np.ones_like(cmb_map, dtype=bool)
    if cmb_map.shape != mask.shape:
        raise ValueError("cmb_map and mask must have same shape")

    theta_bins_rad = np.deg2rad(theta_bins_deg)
    nb = len(theta_bins_deg) - 1
    acc = np.zeros(nb, dtype=float)
    acc2 = np.zeros(nb, dtype=float)
    count = np.zeros(nb, dtype=int)

    n_used = 0
    idxs = np.arange(voids.ra_deg.size)
    if max_voids is not None:
        idxs = idxs[:max_voids]

    for i in idxs:
        ra = np.deg2rad(voids.ra_deg[i])
        dec = np.deg2rad(voids.dec_deg[i])
        v0 = hp.ang2vec(np.pi / 2 - dec, ra)  # healpy uses colatitude, longitude

        # Query a disc out to max radius to limit pixel checks
        pix_disc = hp.query_disc(nside, v0, theta_bins_rad[-1], inclusive=False)

        # Filter by mask
        pix_disc = pix_disc[mask[pix_disc]]
        if pix_disc.size == 0:
            continue

        # Compute angular distances for disc pixels
        vecs = np.array(hp.pix2vec(nside, pix_disc)).T
        ang = np.arccos(np.clip(vecs @ v0, -1.0, 1.0))  # (Npix,)

        # Bin and accumulate
        t = cmb_map[pix_disc] * uK_scale
        for b in range(nb):
            m = (ang >= theta_bins_rad[b]) & (ang < theta_bins_rad[b + 1])
            if np.any(m):
                tb = np.mean(t[m])
                acc[b] += tb
                acc2[b] += tb * tb
                count[b] += 1

        n_used += 1

    prof_mean = np.zeros(nb, dtype=float)
    prof_sem = np.zeros(nb, dtype=float)
    for b in range(nb):
        if count[b] > 0:
            prof_mean[b] = acc[b] / count[b]
            var = max(acc2[b] / count[b] - prof_mean[b] ** 2, 0.0)
            prof_sem[b] = np.sqrt(var / max(count[b], 1))
        else:
            prof_mean[b] = np.nan
            prof_sem[b] = np.nan

    theta_centers = 0.5 * (theta_bins_deg[:-1] + theta_bins_deg[1:])
    return StackResult(theta_bins_deg=theta_centers, prof_mean_uK=prof_mean, prof_sem_uK=prof_sem, n_used=n_used)
