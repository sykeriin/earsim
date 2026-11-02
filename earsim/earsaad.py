"""Ear-SAAD (Geirnaert, Kappel, Kidmose 2025; Zenodo 16536441) channel groups and
approximate electrode geometry on the fsaverage template.

The dataset ships NO electrode coordinates. Scalp channels have standard 10-20 names,
so their positions come from the template montage (exact up to cap fit). The ear
sites do not:
  * around-ear: cEEGrid, 9 (left) / 10 (right) electrodes, numbered 1 = top ... 9/10 = bottom
    (experiment manual). Placed on an arc behind the pinna; arc radius/offset are
    guesses, so E2 reports a sensitivity sweep over them.
  * in-ear: 6 electrodes per ear (A, B, C, E, I, T) inside the ear mould. Their layout
    is not documented and the ear canal is not in the Tier-A head model, so they are
    approximated by a small cluster at the canal entrance. Predictions for them are
    expected to be nearly identical across the six; only ear-averaged / cross-ear
    quantities are meaningful.
"""
from __future__ import annotations

import numpy as np

from .head import load_template
from .layouts import snap_to_scalp, standard

SCALP = ["O1", "Pz", "P3", "P7", "CP1", "CP5", "M1", "Cz", "C3", "T7", "FC1", "FC5", "F3",
         "F7", "Fp1", "O2", "P4", "P8", "CP2", "CP6", "M2", "C4", "T8", "FC2", "FC6", "Fz",
         "F4", "F8", "Fp2"]
IN_EAR_L = ["ELA", "ELB", "ELC", "ELE", "ELI", "ELT"]
IN_EAR_R = ["ERA", "ERB", "ERC", "ERE", "ERI", "ERT"]
CEEGRID_L = [f"cEL{i}" for i in range(1, 10)]
CEEGRID_R = [f"cER{i}" for i in range(1, 11)]
EAR_CHANNELS = IN_EAR_L + CEEGRID_L + IN_EAR_R + CEEGRID_R
GROUPS = {"scalp": SCALP, "in_ear_L": IN_EAR_L, "in_ear_R": IN_EAR_R,
          "around_L": CEEGRID_L, "around_R": CEEGRID_R}


def scalp_positions() -> dict[str, np.ndarray]:
    return standard(SCALP)


def ceegrid_positions(radius: float = 0.032, centre_shift: float = 0.012,
                      sweep=(55.0, 305.0)) -> dict[str, np.ndarray]:
    """cEEGrid arcs behind each pinna (1 = top ... n = bottom), snapped to the BEM scalp.

    The arc lies in the (y, z) plane around a centre `centre_shift` behind the
    pre-auricular point, opening anteriorly where the pinna sits.
    """
    head = load_template()
    out = {}
    for side, key, names in (("L", "lpa", CEEGRID_L), ("R", "rpa", CEEGRID_R)):
        c = head.fiducials[key].copy()
        c[1] -= centre_shift
        ang = np.deg2rad(np.linspace(*sweep, len(names)))
        sgn = -1.0 if side == "L" else 1.0
        pts = np.stack([np.full(len(names), c[0] + sgn * 0.012),
                        c[1] + radius * np.cos(ang),
                        c[2] + radius * np.sin(ang)], axis=1)
        for n, p in zip(names, snap_to_scalp(pts)):
            out[n] = p
    return out


def in_ear_positions(spread: float = 0.006, back: float = 0.010) -> dict[str, np.ndarray]:
    """Cluster of six points at the canal entrance (a ring of radius `spread`)."""
    head = load_template()
    out = {}
    for key, names in (("lpa", IN_EAR_L), ("rpa", IN_EAR_R)):
        c = head.fiducials[key] + np.array([0.0, -back, -0.004])
        ang = np.linspace(0, 2 * np.pi, len(names), endpoint=False)
        pts = np.stack([np.full(len(names), c[0]),
                        c[1] + spread * np.cos(ang), c[2] + spread * np.sin(ang)], axis=1)
        for n, p in zip(names, snap_to_scalp(pts)):
            out[n] = p
    return out


def all_positions(**ceegrid_kw) -> dict[str, np.ndarray]:
    return {**scalp_positions(), **in_ear_positions(), **ceegrid_positions(**ceegrid_kw)}
