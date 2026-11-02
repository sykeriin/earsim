"""E3: zero-shot auditory attention decoding from synthetic ear-EEG (Ear-SAAD preprocessed data).

For every subject s and every held-out trial k:
  real_ear     decoder trained on s's REAL around-ear EEG from the other 5 trials (ceiling)
  synth_krig   decoder trained on SYNTHETIC around-ear EEG made from s's scalp EEG (other trials)
               with the kriging model fitted on the OTHER subjects' scalp covariances
  synth_spline same with spherical-spline interpolation
  scalp        decoder trained and tested on s's scalp EEG (reference: not an ear)
All are tested on s's real around-ear EEG of trial k (scalp for `scalp`). Synthetic methods
never see a real ear sample in training. Channels fed to the decoder: the 19 cEEGrid
electrodes, each re-referenced to the mean of its own ear's grid (reference-free).
Attention is decided per window by comparing the reconstruction's correlation with the
two stream envelopes.
"""
from __future__ import annotations

import argparse
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

from earsim import aad, earsaad, kriging as kr, sources
from earsim.baselines import spherical_spline

warnings.filterwarnings("ignore")
WINS = (10, 30, 60)
AROUND = earsaad.CEEGRID_L + earsaad.CEEGRID_R
EAR = earsaad.EAR_CHANNELS


class Synth:
    """Predict the 31 ear potentials from scalp samples with arbitrary missing channels."""

    def __init__(self, kind, pos, src_names, theta=None, alpha=1e-4, S_all=None):
        self.kind, self.alpha = kind, alpha
        self.src_names = src_names
        self.pos = pos
        self.n = len(src_names)
        self.theta, self.S_all = theta, S_all
        self.cache = {}
        if kind == "krig":
            Sig = np.tensordot(theta[:-1], S_all, 1) + theta[-1] * np.eye(S_all.shape[1])
            self.Sig = Sig

    def W(self, A: np.ndarray):
        key = A.tobytes()
        if key not in self.cache:
            idx = np.where(A)[0]
            if self.kind == "krig":
                self.cache[key] = kr.predict_operator(self.Sig, idx, np.arange(self.n, self.n + len(EAR)))
            else:
                kp = {self.src_names[i]: self.pos[self.src_names[i]] for i in idx}
                tp = {c: self.pos[c] for c in EAR}
                self.cache[key] = spherical_spline(kp, tp, alpha=self.alpha)
        return self.cache[key]

    def __call__(self, X: np.ndarray) -> np.ndarray:
        """X (T, n) scalp with NaN -> (T, 31) ear potentials (avg-ref of available sources)."""
        M = np.isfinite(X)
        Y = np.full((X.shape[0], len(EAR)), np.nan)
        pats, inv = np.unique(M, axis=0, return_inverse=True)
        for pi, A in enumerate(pats):
            if A.sum() < 12:
                continue
            rows = np.where(inv.ravel() == pi)[0]
            Xa = X[np.ix_(rows, np.where(A)[0])]
            Xa = Xa - Xa.mean(1, keepdims=True)
            Y[rows] = Xa @ self.W(A).T
        return Y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default="1-15")
    ap.add_argument("--spline-alpha", type=float, default=1e-4)
    ap.add_argument("--out", default="results/e3")
    ap.add_argument("--ref", choices=["scalp_mean", "local", "cross"], default="cross",
                    help="scalp_mean: ear channel minus the mean of the available scalp channels "
                         "(the shared common reference cancels; NOT device-realistic: the input then contains "
                         "scalp information); local: minus the mean of its own ear grid; cross: minus the "
                         "mean of the opposite ear's grid (device-realistic)")
    args = ap.parse_args()
    global REF
    REF = args.ref
    a, _, b = args.subjects.partition("-")
    subs = list(range(int(a), int(b or a) + 1))

    data = {s: aad.load_prep(s) for s in subs}
    names = data[subs[0]][0]
    nidx = {n: i for i, n in enumerate(names)}
    src_names = earsaad.SCALP
    src_i = [nidx[c] for c in src_names]
    ear_i = [nidx[c] for c in EAR]
    pos = {c: p for c, p in earsaad.all_positions().items() if c in src_names + EAR}
    order = src_names + EAR
    pos = {c: pos[c] for c in order}
    bn, S_raw = sources.cortical_raw(pos)
    S_all = S_raw / sources.avgref_trace(S_raw)[:, None, None]
    n_src = len(src_names)
    B = kr.mean_zero_basis(n_src)
    S_B = np.stack([B.T @ s[:n_src, :n_src] @ B for s in S_all])

    # scalp covariance (1-9 Hz) per subject, trace-normalised
    covs = {}
    for s in subs:
        X = np.vstack([t[:, src_i] for t in data[s][1]])
        P = np.eye(n_src) - 1.0 / n_src
        C = aad.pairwise_cov(X)
        C = P @ C @ P
        covs[s] = C / np.trace(C)
    around_groups = [[AROUND.index(c) for c in earsaad.CEEGRID_L],
                     [AROUND.index(c) for c in earsaad.CEEGRID_R]]
    ear_pos_in_EAR = [EAR.index(c) for c in AROUND]
    real_around_i = [nidx[c] for c in AROUND]

    rows, t0 = [], time.time()
    for s in subs:
        names_s, trials, envs, att = data[s]
        th = kr.fit_theta(S_B, np.stack([B.T @ covs[o] @ B for o in subs if o != s]))
        synth = {"synth_krig": Synth("krig", pos, src_names, theta=th, S_all=S_all),
                 "synth_spline": Synth("spline", pos, src_names, alpha=args.spline_alpha)}
        n_tr = len(trials)
        # precompute per-trial inputs
        scalp = [t[:, src_i] for t in trials]
        if args.ref in ("local", "cross"):
            rf = aad.local_reference if args.ref == "local" else aad.cross_reference
            real_ear = [rf(t[:, real_around_i], around_groups) for t in trials]
            syn_ear = {k: [rf(f(x)[:, ear_pos_in_EAR], around_groups) for x in scalp]
                       for k, f in synth.items()}
        else:
            # real: ear minus mean of available scalp channels; synthetic is already relative to that mean
            real_ear = [t[:, real_around_i] - np.nanmean(t[:, src_i], 1, keepdims=True) for t in trials]
            syn_ear = {k: [f(x)[:, ear_pos_in_EAR] for x in scalp] for k, f in synth.items()}
        ys = [aad.zs(envs[i][:, att[i]]) for i in range(n_tr)]
        for k in range(n_tr):
            tr = [i for i in range(n_tr) if i != k]
            # test features: real ear of trial k, scaled by its own channel std
            Ft_real, _ = aad.prep_features(real_ear[k])
            Ft_scalp, _ = aad.prep_features(scalp[k] - np.nanmean(scalp[k], 1, keepdims=True))
            Ft_syn = {k: aad.prep_features(syn_ear[k][kk])[0] for k, kk in (("synth_krig", k),)}
            plans = {"real_ear": (real_ear, Ft_real),
                     "synth_krig_selftest": (syn_ear["synth_krig"], Ft_syn["synth_krig"]),
                     "synth_krig": (syn_ear["synth_krig"], Ft_real),
                     "synth_spline": (syn_ear["synth_spline"], Ft_real),
                     "scalp": ([x - np.nanmean(x, 1, keepdims=True) for x in scalp], Ft_scalp)}
            for name, (train_sets, Ft) in plans.items():
                Fs = [aad.prep_features(train_sets[i])[0] for i in tr]
                w = aad.train_decoder(Fs, [ys[i] for i in tr])
                recon = Ft @ w
                for ws in WINS:
                    ok, tot = aad.window_accuracy(recon, envs[k], att[k], ws)
                    rows.append({"subject": s, "trial": k, "method": name, "win_s": ws,
                                 "correct": ok, "windows": tot,
                                 "r_att": aad.corr(recon, ys[k])})
        print(f"subject {s} done ({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(Path(args.out).with_suffix(".csv"), index=False)
    summarise(df, Path(args.out))


REF = "cross"


def summarise(df, out):
    g = df.groupby(["subject", "method", "win_s"])[["correct", "windows"]].sum().reset_index()
    g["acc"] = g.correct / g.windows
    L = [f"# E3: zero-shot attention decoding from synthetic ear-EEG (Ear-SAAD), reference = {REF}\n",
         ("**scalp_mean reference is NOT device-realistic**: the decoder input then contains the mean scalp "
          "signal, so this row measures scalp information, not ear information.\n" if REF == "scalp_mean" else ""),
         "Accuracy = fraction of windows where the attended stream is picked; chance = 0.5.\n",
         "Per-subject accuracy, mean over subjects (and number of subjects above chance by "
         "binomial test p<0.05).\n"]
    rows = []
    for (m, w), d in g.groupby(["method", "win_s"]):
        sig = sum(binomtest(int(r.correct), int(r.windows), 0.5, alternative="greater").pvalue < 0.05
                  for r in d.itertuples())
        rows.append({"method": m, "win_s": w, "mean_acc": d.acc.mean(), "median_acc": d.acc.median(),
                     "subjects>chance": f"{sig}/{len(d)}"})
    L.append(pd.DataFrame(rows).round(3).to_markdown(index=False) + "\n")
    L.append("\n## Paired Wilcoxon on per-subject accuracy\n\n| comparison | window | mean diff | p |\n|---|---|---|---|")
    for w in sorted(g.win_s.unique()):
        piv = g[g.win_s == w].pivot(index="subject", columns="method", values="acc")
        for a, b in (("synth_krig", "synth_spline"), ("synth_krig", "real_ear"), ("synth_spline", "real_ear"),
                     ("synth_krig_selftest", "real_ear")):
            d = (piv[a] - piv[b]).dropna()
            p = wilcoxon(d).pvalue if np.any(d != 0) else np.nan
            L.append(f"| {a} - {b} | {w}s | {d.mean():+.3f} | {p:.3g} |")
        d = (piv["synth_krig"] - 0.5)
        L.append(f"| synth_krig - chance | {w}s | {d.mean():+.3f} | {wilcoxon(d).pvalue:.3g} |")
    out.with_suffix(".md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
