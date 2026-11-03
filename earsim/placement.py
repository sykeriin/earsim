"""Electrode-placement search on real ear-EEG for attention decoding.

Candidates are bipolar pairs inside one ear's device (cEEGrid or in-ear mould), which is
what a single earbud / ear-worn strip can record without a distant reference. For every
subject and trial we precompute, once, the Gram matrix G = F^T F and the cross term
b = F^T y of the lagged features of ALL candidate pairs against the attended envelope.
Any subset's subject-specific ridge decoder (trained on the other trials, tested on the
left-out one) is then a small sub-matrix solve, so greedy search is cheap.

Selection is done on the other subjects only (leave-one-subject-out), so the held-out
subject's data never influences which electrodes are chosen.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from . import aad

L = aad.LAGS


def pair_columns(names: list[str], group: list[str]) -> tuple[list[tuple[str, str]], np.ndarray]:
    """All unordered bipolar pairs inside `group` and their column indices in `names`."""
    idx = {n: i for i, n in enumerate(names)}
    pairs = list(combinations(group, 2))
    return pairs, np.array([[idx[a], idx[b]] for a, b in pairs])


@dataclass
class SubjectStats:
    G: list[np.ndarray]      # per trial (P*L, P*L)
    b: list[np.ndarray]      # per trial (P*L,)
    n: list[int]             # samples per trial
    scale: np.ndarray        # (P,) per-pair std used to normalise
    P: int


def pair_features(trial: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """(T, P) bipolar signals a - b, NaN where either channel is missing."""
    return trial[:, cols[:, 0]] - trial[:, cols[:, 1]]


def subject_stats(trials, envs, att, cols: np.ndarray) -> SubjectStats:
    P = len(cols)
    feats = [pair_features(t, cols) for t in trials]
    scale = np.nanstd(np.vstack(feats), 0)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    G, b, n = [], [], []
    for f, e, a in zip(feats, envs, att):
        F = aad.lagged(np.nan_to_num(f / scale))
        y = aad.zs(e[:, a])
        G.append(F.T @ F)
        b.append(F.T @ y)
        n.append(len(y))
    return SubjectStats(G, b, n, scale, P)


def subset_columns(subset: tuple[int, ...], P: int) -> np.ndarray:
    return np.array([k * P + j for k in range(L) for j in subset])


def loto_score(st: SubjectStats, subset: tuple[int, ...], c: float = 0.1) -> float:
    """Mean over left-out trials of the uncentred correlation (cosine similarity) between the
    reconstruction and the attended envelope for this subset. Features are band-passed and
    near zero-mean, so this tracks Pearson r closely; it is only used to rank subsets (final
    window accuracy uses Pearson via aad.corr)."""
    cols = subset_columns(subset, st.P)
    Gs = [g[np.ix_(cols, cols)] for g in st.G]
    bs = [b[cols] for b in st.b]
    tot_G, tot_b = sum(Gs), sum(bs)
    out = []
    for k in range(len(Gs)):
        Gtr, btr = tot_G - Gs[k], tot_b - bs[k]
        if np.trace(Gtr) < 1e-9 or np.trace(Gs[k]) < 1e-9:   # dead channel(s): no usable signal
            out.append(0.0)
            continue
        lam = c * np.trace(Gtr) / Gtr.shape[0]
        w = np.linalg.solve(Gtr + lam * np.eye(len(cols)), btr)
        var = w @ Gs[k] @ w
        out.append(float(bs[k] @ w / np.sqrt(max(var, 1e-12) * st.n[k])))
    return float(np.mean(out))


class ScoreCache:
    """Per-subset, per-subject LOTO scores, shared across all held-out subjects."""

    def __init__(self, stats: dict[int, SubjectStats], c: float = 0.1):
        self.stats, self.c, self.cache = stats, c, {}

    def __call__(self, subset: tuple[int, ...]) -> dict[int, float]:
        key = tuple(sorted(subset))
        if key not in self.cache:
            self.cache[key] = {s: loto_score(st, key, self.c) for s, st in self.stats.items()}
        return self.cache[key]

    def mean_excluding(self, subset, exclude: int) -> float:
        d = self(subset)
        return float(np.mean([v for s, v in d.items() if s != exclude]))


def greedy_select(score: ScoreCache, n_pairs: int, k: int, exclude: int) -> list[int]:
    """Forward selection of k pairs maximising the mean LOTO correlation over subjects != exclude."""
    chosen: list[int] = []
    for _ in range(k):
        best = max((j for j in range(n_pairs) if j not in chosen),
                   key=lambda j: score.mean_excluding(tuple(chosen + [j]), exclude))
        chosen.append(best)
    return chosen


def window_accuracy_for_subset(trials, envs, att, cols, scale, subset, win_s, c=0.1) -> tuple[int, int]:
    """Subject-specific decoder on `subset`, leave-one-trial-out window accuracy."""
    P = len(cols)
    sub = np.array(subset)
    feats = [aad.lagged(np.nan_to_num(pair_features(t, cols)[:, sub] / scale[sub])) for t in trials]
    ys = [aad.zs(e[:, a]) for e, a in zip(envs, att)]
    ok = tot = 0
    for k in range(len(trials)):
        tr = [i for i in range(len(trials)) if i != k]
        F, y = np.vstack([feats[i] for i in tr]), np.concatenate([ys[i] for i in tr])
        if not np.any(F) or not np.any(feats[k]):             # dead channel(s): chance level
            n = sum(1 for _ in range(0, len(envs[k]) - int(win_s * aad.FS) + 1, int(win_s * aad.FS)))
            ok, tot = ok + 0.5 * n, tot + n
            continue
        w = aad.fit_ridge(F, y, c)
        o, t = aad.window_accuracy(feats[k] @ w, envs[k], att[k], win_s)
        ok, tot = ok + o, tot + t
    return ok, tot
