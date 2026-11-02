"""Non-physics baselines that the source-based transfer has to beat."""
from __future__ import annotations

import numpy as np
from mne.channels.interpolation import _make_interpolation_matrix


def spherical_spline(src_pos: dict[str, np.ndarray], tgt_pos: dict[str, np.ndarray],
                     alpha: float | None = 1e-5) -> np.ndarray:
    """Perrin et al. (1989) spherical-spline interpolation matrix (n_tgt x n_src).

    Positions are centred on a least-squares sphere fitted to the source electrodes
    and normalised to it. `alpha` is the diagonal regulariser; MNE's default 1e-5
    smooths heavily (it does not reproduce the input channels), so treat it as a
    hyperparameter and tune it on dev data like the inverse SNR.
    """
    S = np.array(list(src_pos.values()))
    D = np.array(list(tgt_pos.values()))
    centre = _fit_sphere_centre(S)
    return _make_interpolation_matrix(S - centre, D - centre, alpha=alpha)


def nearest_neighbour(src_pos, tgt_pos) -> np.ndarray:
    S = np.array(list(src_pos.values()))
    D = np.array(list(tgt_pos.values()))
    W = np.zeros((len(D), len(S)))
    W[np.arange(len(D)), np.argmin(np.linalg.norm(D[:, None] - S[None], axis=2), axis=1)] = 1
    return W


def _fit_sphere_centre(P: np.ndarray) -> np.ndarray:
    # linear least squares: |p|^2 = 2 p.c + (r^2 - |c|^2)
    A = np.c_[2 * P, np.ones(len(P))]
    b = (P**2).sum(1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    return sol[:3]
