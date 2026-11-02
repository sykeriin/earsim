"""Electrode layouts in fsaverage MRI coordinates (metres).

Wearable layouts here are geometric approximations defined relative to the
template's pre-auricular points (LPA/RPA) and snapped to the BEM scalp. They are
not measured device geometries. Replace with measured coordinates (e.g. from a
CAD file or photogrammetry) before drawing device-specific conclusions.
"""
from __future__ import annotations

import mne
import numpy as np

from .head import load_template

# MRI frame of fsaverage: +x right, +y anterior, +z superior.


def standard(names: list[str] | None = None, kind: str = "standard_1005") -> dict[str, np.ndarray]:
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        pos = mne.channels.make_standard_montage(kind).get_positions()
    assert pos["coord_frame"] == "mri"
    ch = pos["ch_pos"]
    if names is None:
        return dict(ch)
    return {n: np.asarray(ch[n]) for n in names}


def snap_to_scalp(points: np.ndarray) -> np.ndarray:
    head = load_template()
    rr = head.scalp_rr
    idx = np.argmin(np.linalg.norm(rr[None] - points[:, None], axis=2), axis=1)
    return rr[idx]


def _side_sign(side: str) -> float:
    return -1.0 if side == "L" else 1.0


def ceegrid_like(side: str = "L", n: int = 10, radius: float = 0.030,
                 posterior_shift: float = 0.012) -> dict[str, np.ndarray]:
    """C-shaped array behind the ear, cEEGrid-like (10 electrodes per side).

    Arc in the (y, z) plane centred slightly behind the pre-auricular point,
    opening anteriorly (the gap is where the pinna attaches), then snapped to scalp.
    """
    head = load_template()
    c = head.fiducials["lpa" if side == "L" else "rpa"].copy()
    c[1] -= posterior_shift
    # angles measured from +y (anterior); open arc from 60 deg to 300 deg
    ang = np.deg2rad(np.linspace(60, 300, n))
    pts = np.stack([np.full(n, c[0] + _side_sign(side) * 0.01),
                    c[1] + radius * np.cos(ang),
                    c[2] + radius * np.sin(ang)], axis=1)
    pts = snap_to_scalp(pts)
    return {f"cEG{side}{i + 1:02d}": p for i, p in enumerate(pts)}


def in_ear_like(side: str = "L") -> dict[str, np.ndarray]:
    """Crude in-ear proxies: canal entrance and concha.

    Tier A cannot model the ear canal; these points sit on the BEM scalp right at
    the pre-auricular region, i.e. they are an "around the tragus" approximation.
    """
    head = load_template()
    c = head.fiducials["lpa" if side == "L" else "rpa"]
    offs = {"canal": np.array([0.0, -0.008, -0.004]),
            "concha": np.array([0.0, -0.015, 0.006])}
    pts = snap_to_scalp(np.stack([c + o for o in offs.values()]))
    return {f"IE{side}_{k}": p for k, p in zip(offs, pts)}


def ear_candidates(radius: float = 0.040, min_spacing: float = 0.008) -> dict[str, np.ndarray]:
    """Candidate electrode sites on the scalp within `radius` of each ear.

    Greedy farthest-point thinning to `min_spacing`. Used by the placement optimiser.
    """
    head = load_template()
    out = {}
    for side, key in (("L", "lpa"), ("R", "rpa")):
        d = np.linalg.norm(head.scalp_rr - head.fiducials[key], axis=1)
        pts = head.scalp_rr[d < radius]
        keep = [pts[np.argmin(np.linalg.norm(pts - head.fiducials[key], axis=1))]]
        for p in pts[np.argsort(np.linalg.norm(pts - head.fiducials[key], axis=1))]:
            if min(np.linalg.norm(np.array(keep) - p, axis=1)) >= min_spacing:
                keep.append(p)
        out.update({f"C{side}{i:03d}": p for i, p in enumerate(keep)})
    return out


def wearable_layouts() -> dict[str, dict[str, np.ndarray]]:
    """Named target layouts used in experiments."""
    lay = {}
    lay["ceegrid"] = {**ceegrid_like("L"), **ceegrid_like("R")}
    lay["in_ear"] = {**in_ear_like("L"), **in_ear_like("R")}
    # headband: ring at roughly the FT/T/TP line plus forehead (Fp) sites
    lay["headband"] = standard(["Fp1", "Fpz", "Fp2", "AF7", "AF8", "FT9", "FT10", "T9", "T10"])
    return lay
