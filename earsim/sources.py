"""Augmented source space for the covariance model.

Every source group k contributes one basis covariance S_k = A_k A_k^T over the
full channel set, where A_k holds the lead-field columns of the group:

  * cortex: one group per aparc parcel (vertices within a parcel share one variance),
  * ocular: corneo-retinal dipoles, common-mode (L+R) and differential (L-R) patterns,
  * muscle: temporalis and masseter dipoles on each side.

Cortical columns come from the BEM forward. Ocular and muscle sources lie OUTSIDE
the inner skull, where MNE's BEM cannot place sources, so they use an analytic
infinite-homogeneous-medium dipole field. That ignores the skull and the head
boundary. It is a deliberate approximation: only the spatial pattern matters
(amplitudes are fitted), but near the source the true pattern differs. It is the
main modelling risk of the approach and is ablated against cortex-only.

All S_k are scaled to unit trace after average referencing, so fitted variances are
comparable across groups.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import mne
import numpy as np

from .head import SUBJECTS_DIR, fixed, forward, load_template, make_info

# eyeball centre relative to the nasion (m), MRI frame: lateral, anterior, vertical
EYE_OFFSET = np.array([0.032, -0.004, -0.006])


def _unit_trace(S: np.ndarray) -> np.ndarray:
    n = S.shape[0]
    P = np.eye(n) - 1.0 / n
    t = np.trace(P @ S @ P)
    return S / t


def cortical_bases(ch_pos: dict[str, np.ndarray], parc: str = "aparc"):
    """(names, S[K, n, n]) for one group per parcel, from the BEM lead field."""
    info = make_info(ch_pos)
    fwd = fixed(forward(info))
    G = fwd["sol"]["data"]  # (n_ch, n_src), average-referenced by MNE
    vert = [s["vertno"] for s in fwd["src"]]
    offsets = np.cumsum([0] + [len(v) for v in vert])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        labels = mne.read_labels_from_annot("fsaverage", parc, subjects_dir=SUBJECTS_DIR,
                                            verbose=False)
    names, mats = [], []
    for lab in labels:
        if lab.name.startswith(("unknown", "corpuscallosum")):
            continue
        hemi = 0 if lab.hemi == "lh" else 1
        idx = np.nonzero(np.isin(vert[hemi], lab.vertices))[0] + offsets[hemi]
        if len(idx) == 0:
            continue
        A = G[:, idx]
        names.append(f"ctx:{lab.name}")
        mats.append(_unit_trace(A @ A.T))
    return names, np.stack(mats)


def dipole_field(ch_xyz: np.ndarray, r0: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Infinite-medium potential of a unit dipole p at r0, evaluated at ch_xyz (n,3)."""
    d = ch_xyz - r0
    return d @ p / (4 * np.pi * np.linalg.norm(d, axis=1) ** 3)


def _inward(point: np.ndarray, depth: float) -> np.ndarray:
    """Move a scalp point `depth` metres towards the head centre."""
    head = load_template()
    centre = head.scalp_rr.mean(0)
    v = centre - point
    return point + depth * v / np.linalg.norm(v)


def extra_source_positions() -> dict[str, np.ndarray]:
    """Template positions (MRI frame, metres) of ocular and muscle generators."""
    head = load_template()
    nas = head.fiducials["nasion"]
    pos = {"eye_L": nas + EYE_OFFSET * np.array([-1, 1, 1]),
           "eye_R": nas + EYE_OFFSET}
    rr = head.scalp_rr

    def snap(p):
        return rr[np.argmin(np.linalg.norm(rr - p, axis=1))]

    for side, key, sgn in (("L", "lpa", -1.0), ("R", "rpa", 1.0)):
        e = head.fiducials[key]
        pos[f"temporalis_{side}"] = _inward(snap(e + np.array([0.0, 0.020, 0.040])), 0.006)
        pos[f"masseter_{side}"] = _inward(snap(e + np.array([0.0, 0.035, -0.035])), 0.006)
    return pos


def extra_bases(ch_pos: dict[str, np.ndarray]):
    """(names, S[K, n, n]) for ocular (6) and muscle (12) groups."""
    X = np.array(list(ch_pos.values()))
    pos = extra_source_positions()
    axes = {"x": np.array([1.0, 0, 0]), "y": np.array([0, 1.0, 0]), "z": np.array([0, 0, 1.0])}
    names, mats = [], []
    gL = {a: dipole_field(X, pos["eye_L"], v) for a, v in axes.items()}
    gR = {a: dipole_field(X, pos["eye_R"], v) for a, v in axes.items()}
    for a in axes:
        for tag, g in (("common", gL[a] + gR[a]), ("diff", gL[a] - gR[a])):
            names.append(f"eye:{tag}_{a}")
            mats.append(_unit_trace(np.outer(g, g)))
    for name in ("temporalis_L", "temporalis_R", "masseter_L", "masseter_R"):
        for a, v in axes.items():
            g = dipole_field(X, pos[name], v)
            names.append(f"muscle:{name}_{a}")
            mats.append(_unit_trace(np.outer(g, g)))
    return names, np.stack(mats)


def all_bases(ch_pos: dict[str, np.ndarray], use_extra: bool = True):
    """Basis covariances over the channels in `ch_pos` (reference-free)."""
    names, S = cortical_bases(ch_pos)
    if use_extra:
        n2, S2 = extra_bases(ch_pos)
        names, S = names + n2, np.concatenate([S, S2])
    return names, S
