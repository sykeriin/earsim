"""Source-space transfer: scalp recording -> sources -> arbitrary electrode layout.

For any *linear* inverse (MNE, eLORETA) the whole pipeline collapses to one
matrix T = R_tgt @ G_tgt @ K, where K is the inverse kernel for the recorded
montage and R_tgt is the target reference operator. T has rank <= n_scalp, so
synthesis is physics-informed spatial extrapolation, not new information. Sparse
(non-linear) inverses are handled per window in `sparse.py` (not yet built).

dSPM and sLORETA are deliberately excluded: they output noise-normalised
statistics, not current estimates, so re-projecting them through a lead field
is not physically meaningful.
"""
from __future__ import annotations

from dataclasses import dataclass

import mne
import numpy as np

from . import cache
from .head import fixed, forward, leadfield, make_info

LINEAR_METHODS = ("MNE", "eLORETA")


def avg_ref_operator(n: int) -> np.ndarray:
    return np.eye(n) - np.ones((n, n)) / n


@dataclass
class Transfer:
    T: np.ndarray            # (n_target, n_source_channels)
    src_names: list[str]
    tgt_names: list[str]
    method: str
    lambda2: float

    def apply(self, data: np.ndarray) -> np.ndarray:
        """data: (n_source_channels, n_times), any reference. Avg-ref is applied inside K."""
        return self.T @ data


def inverse_kernel(src_pos: dict[str, np.ndarray], method: str = "MNE", snr: float = 3.0,
                   depth: float | None = 0.8, sfreq: float = 160.0):
    """Explicit inverse kernel K (n_sources x n_src_channels), fixed orientation.

    Includes the average-reference projection and whitening with an ad-hoc
    (diagonal) noise covariance, so K expects data in any common reference.
    """
    if method not in LINEAR_METHODS:
        raise ValueError(f"method must be one of {LINEAR_METHODS}, got {method}")
    info0 = make_info(src_pos, sfreq)
    fwd_free = forward(info0)
    fwd = fixed(fwd_free)
    n = len(src_pos)
    ev = mne.EvokedArray(np.eye(n), info0, tmin=0.0, verbose=False)
    ev.set_eeg_reference("average", projection=True, verbose=False)
    info = ev.info
    cov = mne.make_ad_hoc_cov(info, verbose=False)
    inv = mne.minimum_norm.make_inverse_operator(info, fwd_free, cov, loose=0.0, depth=depth,
                                                 fixed=True, verbose=False)
    lambda2 = 1.0 / snr**2
    stc = mne.minimum_norm.apply_inverse(ev, inv, lambda2=lambda2, method=method,
                                         pick_ori=None, verbose=False)
    return stc.data, fwd, lambda2


def build_transfer(src_pos: dict[str, np.ndarray], tgt_pos: dict[str, np.ndarray],
                   method: str = "MNE", snr: float = 3.0, depth: float | None = 0.8,
                   tgt_reference: str | np.ndarray = "average_of_source") -> Transfer:
    """Build T mapping recorded channels -> target electrodes.

    tgt_reference:
      "average_of_source": targets referenced to the mean of the *source* channels'
                           modelled potentials (matches avg-ref'd recorded data; used
                           for leave-channels-out validation).
      "infinity":          reference-free potentials.
      "average":           common average over target channels.
      ndarray (n_out, n_tgt): any explicit reference/montage operator, e.g. bipolar.
    """
    ref_key = tgt_reference if isinstance(tgt_reference, str) else \
        cache.key(None, ref=np.asarray(tgt_reference).round(9).tolist())
    k = cache.key(src_pos, tgt=cache.key(tgt_pos), method=method, snr=snr, depth=depth,
                  ref=ref_key, v=1)
    f = cache.path("T", k, "npz")
    if f.exists():
        z = np.load(f, allow_pickle=False)
        return Transfer(z["T"], list(src_pos), list(z["names"]), method, float(z["lambda2"]))
    tr = _build_transfer(src_pos, tgt_pos, method, snr, depth, tgt_reference)
    np.savez(f, T=tr.T, names=np.array(tr.tgt_names), lambda2=tr.lambda2)
    return tr


def _build_transfer(src_pos, tgt_pos, method, snr, depth, tgt_reference) -> Transfer:
    K, fwd_src, lambda2 = inverse_kernel(src_pos, method, snr, depth)
    G_tgt, fwd_tgt = leadfield(make_info(tgt_pos))
    _check_same_sources(fwd_src, fwd_tgt)
    X = G_tgt @ K  # reference-free target potentials per unit source-channel input
    names = list(tgt_pos)
    if isinstance(tgt_reference, np.ndarray):
        X = tgt_reference @ X
        names = [f"ref{i}" for i in range(X.shape[0])]
    elif tgt_reference == "average_of_source":
        G_src = fwd_src["sol"]["data"]
        X = X - (G_src @ K).mean(axis=0, keepdims=True)
    elif tgt_reference == "average":
        X = avg_ref_operator(X.shape[0]) @ X
    elif tgt_reference != "infinity":
        raise ValueError(tgt_reference)
    return Transfer(X, list(src_pos), names, method, lambda2)


def _check_same_sources(a: mne.Forward, b: mne.Forward) -> None:
    va = [s["vertno"] for s in a["src"]]
    vb = [s["vertno"] for s in b["src"]]
    if any(len(x) != len(y) or not np.array_equal(x, y) for x, y in zip(va, vb)):
        raise RuntimeError("source spaces differ between forwards (mindist pruning?)")


def bipolar(names: list[str], pairs: list[tuple[str, str]]) -> np.ndarray:
    """Reference operator for bipolar derivations a - b."""
    idx = {n: i for i, n in enumerate(names)}
    R = np.zeros((len(pairs), len(names)))
    for k, (a, b) in enumerate(pairs):
        R[k, idx[a]] = 1.0
        R[k, idx[b]] = -1.0
    return R
