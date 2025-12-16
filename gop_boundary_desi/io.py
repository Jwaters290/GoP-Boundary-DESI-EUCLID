from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from astropy.io import fits

@dataclass(frozen=True)
class VoidCatalog:
    ra_deg: np.ndarray
    dec_deg: np.ndarray
    z: np.ndarray
    r_mpc_h: Optional[np.ndarray] = None

def _require_cols(df: pd.DataFrame, cols: Tuple[str, ...]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

def read_void_catalog(path: str | Path) -> VoidCatalog:
    """
    Read a void catalog from CSV/Parquet or FITS.

    Required columns: RA (deg), DEC (deg), Z.
    Optional: R (void radius, e.g., Mpc/h).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in [".csv"]:
        df = pd.read_csv(path)
    elif path.suffix.lower() in [".parquet"]:
        df = pd.read_parquet(path)
    elif path.suffix.lower() in [".fits", ".fit", ".fz"]:
        with fits.open(path) as hdul:
            # Try first table HDU with data
            for hdu in hdul:
                if hasattr(hdu, "data") and hdu.data is not None and len(hdu.data) > 0:
                    df = pd.DataFrame(np.array(hdu.data).byteswap().newbyteorder())
                    break
            else:
                raise ValueError(f"No table data found in FITS: {path}")
    else:
        raise ValueError(f"Unsupported file type: {path.suffix}")

    # Normalize column names (common variants)
    colmap = {c: c.strip().upper() for c in df.columns}
    df = df.rename(columns=colmap)

    # Common aliases
    aliases = {
        "RA": ("RA", "RA_DEG", "ALPHA", "ALPHA_J2000"),
        "DEC": ("DEC", "DEC_DEG", "DELTA", "DELTA_J2000"),
        "Z": ("Z", "REDSHIFT"),
        "R": ("R", "R_VOID", "RADIUS", "R_MPC_H", "R_MPC")
    }

    def pick(name: str) -> Optional[str]:
        for a in aliases[name]:
            if a in df.columns:
                return a
        return None

    ra_col = pick("RA")
    dec_col = pick("DEC")
    z_col = pick("Z")
    if not (ra_col and dec_col and z_col):
        raise ValueError(f"Could not find RA/DEC/Z in columns: {list(df.columns)}")

    r_col = pick("R")

    ra = df[ra_col].to_numpy(dtype=float)
    dec = df[dec_col].to_numpy(dtype=float)
    z = df[z_col].to_numpy(dtype=float)
    r = df[r_col].to_numpy(dtype=float) if r_col else None

    mask = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
    if r is not None:
        mask &= np.isfinite(r)

    ra, dec, z = ra[mask], dec[mask], z[mask]
    r = r[mask] if r is not None else None
    return VoidCatalog(ra_deg=ra, dec_deg=dec, z=z, r_mpc_h=r)
