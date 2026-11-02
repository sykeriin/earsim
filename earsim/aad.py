"""Auditory attention decoding helpers for the Ear-SAAD preprocessed files (20 Hz, 1-9 Hz).

Backward model: ridge regression from time-lagged EEG (0-250 ms ahead) to the attended
speech envelope; a window is classified by which of the two stream envelopes correlates
better with the reconstruction.
"""
from __future__ import annotations

import numpy as np
import scipy.io as sio

from .head import DATA_DIR

PREP = DATA_DIR / "earsaad/preprocessedData"
FS = 20
LAGS = 6  # 0..5 samples = 0..250 ms


def load_prep(n: int):
    m = sio.loadmat(PREP / f"dataSubject{n}.mat")
    names = [str(c[0]) for c in m["channels"].ravel()]
    trials = [np.asarray(m["eegTrials"][i, 0], float) for i in range(m["eegTrials"].shape[0])]
    env = [np.asarray(m["envelopes"][i, 0], float) for i in range(m["envelopes"].shape[0])]
    att = m["attSpeaker"].ravel().astype(int) - 1  # 0 or 1: which stream is attended
    return names, trials, env, att


def pairwise_cov(X: np.ndarray) -> np.ndarray:
    """Covariance of (T, n) data with NaNs by pairwise-complete estimation, made PSD."""
    M = np.isfinite(X).astype(float)
    Z = np.nan_to_num(X - np.nanmean(X, 0))
    C = (Z.T @ Z) / np.maximum(M.T @ M, 1)
    w, V = np.linalg.eigh((C + C.T) / 2)
    return (V * np.maximum(w, 1e-8 * w.max())) @ V.T


def local_reference(Y: np.ndarray, groups: list[list[int]]) -> np.ndarray:
    """Subtract each group's per-sample mean over its available (finite) channels."""
    out = np.full((Y.shape[0], sum(len(g) for g in groups)), np.nan)
    col = 0
    for g in groups:
        sub = Y[:, g]
        n_ok = np.isfinite(sub).sum(1)
        mean = np.where(n_ok > 2, np.nanmean(np.where(np.isfinite(sub), sub, np.nan), 1), np.nan)
        out[:, col:col + len(g)] = sub - mean[:, None]
        col += len(g)
    return out


def cross_reference(Y: np.ndarray, groups: list[list[int]]) -> np.ndarray:
    """Each channel minus the per-sample mean of the OTHER ear's group (device-realistic)."""
    g0, g1 = groups
    out = np.full((Y.shape[0], len(g0) + len(g1)), np.nan)
    for tgt, other, off in ((g0, g1, 0), (g1, g0, len(g0))):
        sub = Y[:, other]
        n_ok = np.isfinite(sub).sum(1)
        mean = np.where(n_ok > 2, np.nanmean(np.where(np.isfinite(sub), sub, np.nan), 1), np.nan)
        out[:, off:off + len(tgt)] = Y[:, tgt] - mean[:, None]
    return out


def lagged(X: np.ndarray) -> np.ndarray:
    """(T, c) -> (T, c*LAGS) with columns [x(t), x(t+1), ..., x(t+LAGS-1)]; zero-padded."""
    T, c = X.shape
    out = np.zeros((T, c * LAGS))
    for k in range(LAGS):
        out[: T - k, k * c:(k + 1) * c] = X[k:]
    return out


def prep_features(X: np.ndarray, scale: np.ndarray | None = None):
    """Zero-fill NaN, scale channels by their std (computed here unless given)."""
    if scale is None:
        scale = np.nanstd(X, 0)
        scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    return lagged(np.nan_to_num(X / scale)), scale


def fit_ridge(F: np.ndarray, y: np.ndarray, c: float) -> np.ndarray:
    G = F.T @ F
    lam = c * np.trace(G) / G.shape[0]
    return np.linalg.solve(G + lam * np.eye(G.shape[0]), F.T @ y)


def zs(x):
    return (x - x.mean(0)) / (x.std(0) + 1e-12)


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


def train_decoder(Fs: list[np.ndarray], ys: list[np.ndarray],
                  grid=(1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)) -> np.ndarray:
    """Ridge decoder with the regularisation chosen by leave-one-trial-out on the training data."""
    best, best_c = -np.inf, grid[0]
    for c in grid:
        sc = []
        for k in range(len(Fs)):
            tr = [i for i in range(len(Fs)) if i != k]
            w = fit_ridge(np.vstack([Fs[i] for i in tr]), np.concatenate([ys[i] for i in tr]), c)
            sc.append(corr(Fs[k] @ w, ys[k]))
        if np.mean(sc) > best:
            best, best_c = np.mean(sc), c
    return fit_ridge(np.vstack(Fs), np.concatenate(ys), best_c)


def window_accuracy(recon: np.ndarray, env: np.ndarray, att: int, win_s: float) -> tuple[int, int]:
    """(# correct, # windows) classifying the attended stream in non-overlapping windows."""
    n = int(win_s * FS)
    ok = tot = 0
    for s in range(0, len(recon) - n + 1, n):
        r = recon[s:s + n]
        c = [corr(r, env[s:s + n, j]) for j in (0, 1)]
        ok += int(np.argmax(c) == att)
        tot += 1
    return ok, tot
