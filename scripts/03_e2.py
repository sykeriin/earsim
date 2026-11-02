"""E2: predict real ear-EEG (Ear-SAAD) from the scalp channels, leave-one-subject-out.

Sources: scalp channels that are good in (almost) all subjects.
Targets: 31 ear channels (12 in-ear, 19 cEEGrid). The dataset has no electrode
coordinates, so ear positions are template approximations (earsim/earsaad.py) and the
result is reported with a geometry sensitivity sweep.

Outputs are reference-free (rows sum to zero) because the common reference is not
known exactly (manual: Fpz/Fp1):
  around_local  each cEEGrid electrode minus the mean of its own ear's cEEGrid
  around_cross  mean(cEEGrid L) - mean(cEEGrid R)
  in_cross      mean(in-ear L) - mean(in-ear R)
  around_in     mean(cEEGrid) - mean(in-ear), per ear (2 outputs)
  in_local      each in-ear electrode minus its own ear's in-ear mean (expected ~unpredictable)

Methods (everything fitted only on the other subjects): spline (alpha tuned), MNE
(identity source prior, SNR tuned), krig_cortex (parcel variances fitted on the other
subjects' scalp covariances), and `oracle_lsq`: a ridge regression scalp->ear learned on
the other subjects' *paired* ear data, a ceiling the physics model cannot be expected to
beat because it never sees ear data.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from earsim import cache, earsaad, kriging as kr, sources
from earsim.baselines import spherical_spline
from earsim.transfer import build_transfer

warnings.filterwarnings("ignore")
SPLINE_ALPHAS = [1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
MNE_SNRS = [1.0, 3.0, 10.0]
EAR = earsaad.EAR_CHANNELS
OUTPUT_SETS = ["around_local", "around_cross", "in_cross", "around_in", "in_local"]


def load_covs():
    files = sorted(Path(cache.CACHE_DIR).glob("cov-earsaad-sub*.npz"))
    subs, C, good, names, bands = [], [], [], None, None
    for f in files:
        z = np.load(f, allow_pickle=False)
        subs.append(int(f.stem[-2:]))
        C.append(z["C"])
        good.append(z["good"])
        names = [str(n) for n in z["ch_names"]]
        bands = [str(b) for b in z["bands"]]
    return np.array(subs), np.stack(C), np.stack(good), names, bands


def output_operator(ear_good: np.ndarray) -> dict[str, np.ndarray]:
    """H[set] (n_out x 31) over the EAR channel order, zero columns for bad channels."""
    idx = {n: i for i, n in enumerate(EAR)}
    g = {n: bool(ear_good[idx[n]]) for n in EAR}

    def avg(chs):
        v = np.zeros(len(EAR))
        ok = [c for c in chs if g[c]]
        if ok:
            v[[idx[c] for c in ok]] = 1.0 / len(ok)
        return v, len(ok)

    H = {}
    rows = []
    for grp in (earsaad.CEEGRID_L, earsaad.CEEGRID_R):
        m, n = avg(grp)
        for c in grp:
            if g[c] and n > 2:
                e = np.zeros(len(EAR)); e[idx[c]] = 1
                rows.append(e - m)
    H["around_local"] = np.array(rows)
    rows = []
    for grp in (earsaad.IN_EAR_L, earsaad.IN_EAR_R):
        m, n = avg(grp)
        for c in grp:
            if g[c] and n > 2:
                e = np.zeros(len(EAR)); e[idx[c]] = 1
                rows.append(e - m)
    H["in_local"] = np.array(rows)
    aL, nAL = avg(earsaad.CEEGRID_L); aR, nAR = avg(earsaad.CEEGRID_R)
    iL, nIL = avg(earsaad.IN_EAR_L); iR, nIR = avg(earsaad.IN_EAR_R)
    H["around_cross"] = (aL - aR)[None] if nAL and nAR else np.zeros((0, len(EAR)))
    H["in_cross"] = (iL - iR)[None] if nIL and nIR else np.zeros((0, len(EAR)))
    rows = []
    if nAL and nIL: rows.append(aL - iL)
    if nAR and nIR: rows.append(aR - iR)
    H["around_in"] = np.array(rows) if rows else np.zeros((0, len(EAR)))
    return H


def score(W, Cs, src_idx, ear_idx, H):
    """Per output-set (mean r, mean gain) of pred = H W P x vs true = H y, from covariance Cs (64x64)."""
    n = len(src_idx)
    P = np.eye(n) - 1.0 / n
    Cx = np.nan_to_num(Cs[np.ix_(src_idx, src_idx)])
    Cxy = np.nan_to_num(Cs[np.ix_(src_idx, ear_idx)])
    Cy = np.nan_to_num(Cs[np.ix_(ear_idx, ear_idx)])
    out = {}
    for k, Hk in H.items():
        if len(Hk) == 0:
            out[k] = (np.nan, np.nan)
            continue
        D = Hk @ W @ P
        cov = np.einsum("oi,io->o", D, Cxy @ Hk.T)
        vp = np.einsum("oi,ij,oj->o", D, Cx, D)
        vt = np.einsum("oi,ij,oj->o", Hk, Cy, Hk)
        ok = (vp > 0) & (vt > 0)
        r = cov[ok] / np.sqrt(vp[ok] * vt[ok])
        out[k] = (float(r.mean()) if ok.any() else np.nan,
                  float(np.sqrt(vp[ok] / vt[ok]).mean()) if ok.any() else np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e2")
    ap.add_argument("--min-subjects", type=int, default=4)
    ap.add_argument("--src-frac", type=float, default=0.7,
                    help="scalp channel is a source if good in at least this fraction of subjects; "
                         "subjects missing any source channel are excluded")
    ap.add_argument("--sweep", action="store_true", help="geometry sensitivity sweep (krig/spline only)")
    ap.add_argument("--bands", default="delta,theta,alpha,beta,gamma,broad")
    args = ap.parse_args()

    subs, C, good, names, all_bands = load_covs()
    bands = args.bands.split(",")
    print(f"{len(subs)} subjects: {subs.tolist()}", flush=True)
    if len(subs) < args.min_subjects:
        raise SystemExit("not enough subjects cached yet")
    nidx = {n: i for i, n in enumerate(names)}
    scalp_i = np.array([nidx[c] for c in earsaad.SCALP])
    frac_good = good[:, scalp_i].mean(0)
    src_names = [c for c, f in zip(earsaad.SCALP, frac_good) if f >= args.src_frac]
    use = np.array([all(good[s, nidx[c]] for c in src_names) for s in range(len(subs))])
    print(f"source set ({len(src_names)}): {src_names}\n  usable subjects: {subs[use].tolist()}", flush=True)
    C, good, subs = C[use], good[use], subs[use]
    S_ = len(subs)
    src_idx = np.array([nidx[c] for c in src_names])
    ear_idx = np.array([nidx[c] for c in EAR])
    Hs = [output_operator(good[s][ear_idx]) for s in range(S_)]

    configs = [dict(radius=0.032, centre_shift=0.012)]
    if args.sweep:
        configs = [dict(radius=r, centre_shift=c) for r in (0.025, 0.032, 0.040) for c in (0.006, 0.012, 0.020)]

    rows, t0 = [], time.time()
    for ci, cfg in enumerate(configs):
        full = ci == 0 and not args.sweep or (args.sweep and cfg == dict(radius=0.032, centre_shift=0.012))
        pos = {c: p for c, p in earsaad.all_positions(**cfg).items() if c in src_names + EAR}
        order = src_names + EAR
        pos = {c: pos[c] for c in order}
        kpos = {c: pos[c] for c in src_names}
        tpos = {c: pos[c] for c in EAR}
        n_src = len(src_names)
        _, S_all = sources.cortical_bases(pos)
        S_src = S_all[:, :n_src, :n_src]
        B = kr.mean_zero_basis(n_src)
        S_B = np.stack([B.T @ s @ B for s in S_src])
        ops_static = {f"spline|a={a}": spherical_spline(kpos, tpos, alpha=a) for a in SPLINE_ALPHAS}
        if full:
            for snr in MNE_SNRS:
                ops_static[f"MNE|snr={snr:g}"] = build_transfer(kpos, tpos, method="MNE", snr=snr).T
        print(f"[{ci + 1}/{len(configs)}] {cfg} bases+baselines ({time.time() - t0:.0f}s)", flush=True)

        for bi, band in enumerate(all_bands):
            if band not in bands:
                continue
            Cb = C[:, bi]                                          # (S, 64, 64)
            tr = np.array([np.trace(np.nan_to_num(Cb[s][np.ix_(src_idx, src_idx)])) for s in range(S_)])
            # static operators: score for every subject
            static_scores = {k: [score(W, Cb[s], src_idx, ear_idx, Hs[s]) for s in range(S_)]
                             for k, W in ops_static.items()}
            for t in range(S_):
                train = [s for s in range(S_) if s != t]
                # LOSO selection of alpha / snr on the training subjects
                pick = {}
                for fam in ("spline", "MNE"):
                    keys = [k for k in ops_static if k.startswith(fam)]
                    if not keys:
                        continue
                    def prim(k):
                        vals = [np.nanmean([static_scores[k][s][o][0] for o in OUTPUT_SETS[:4]]) for s in train]
                        return np.nanmean(vals)
                    pick[fam] = max(keys, key=prim)
                # kriging: fit theta on training subjects' scalp covariance
                Ctr = np.stack([B.T @ (np.nan_to_num(Cb[s][np.ix_(src_idx, src_idx)]) / tr[s]) @ B
                                for s in train])
                th = kr.fit_theta(S_B, Ctr)
                Sig = np.tensordot(th[:-1], S_all, 1) + th[-1] * np.eye(S_all.shape[1])
                W_k = kr.predict_operator(Sig, np.arange(n_src), np.arange(n_src, n_src + len(EAR)))
                methods = {"krig_cortex": W_k, **{fam: ops_static[k] for fam, k in pick.items()}}
                if full:
                    # paired-data ceiling: ridge regression on the other subjects' covariances
                    P = np.eye(n_src) - 1.0 / n_src
                    Cxx = sum(P @ np.nan_to_num(Cb[s][np.ix_(src_idx, src_idx)]) @ P / tr[s] for s in train)
                    Cxy = sum(P @ np.nan_to_num(Cb[s][np.ix_(src_idx, ear_idx)]) / tr[s] for s in train)
                    lam = 1e-2 * np.trace(Cxx) / n_src
                    methods["oracle_lsq"] = (np.linalg.solve(Cxx + lam * np.eye(n_src), Cxy)).T
                for name, W in methods.items():
                    sc = score(W, Cb[t], src_idx, ear_idx, Hs[t])
                    for o, (r, g) in sc.items():
                        rows.append({"cfg": f"r={cfg['radius']}|c={cfg['centre_shift']}", "subject": int(subs[t]),
                                     "band": band, "method": name, "outputs": o, "r": r, "gain": g,
                                     "chosen": pick.get(name, "-")})
            print(f"   band {band} done ({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    df.to_csv(out.with_suffix(".csv"), index=False)
    summarise(df, out, subs)


def summarise(df, out, subs):
    L = [f"# E2: scalp -> real ear-EEG (Ear-SAAD), leave-one-subject-out, n={len(subs)} subjects\n",
         "Ear electrode positions are template approximations (no coordinates ship with the dataset).\n"
         "Median over subjects of mean r per output set (reference-free outputs).\n"]
    for cfg in df.cfg.unique():
        d = df[df.cfg == cfg]
        L.append(f"\n## geometry {cfg}\n")
        for o in OUTPUT_SETS:
            t = d[d.outputs == o].groupby(["band", "method"]).r.median().unstack()
            L.append(f"\n**{o}**\n\n" + t.round(3).to_markdown() + "\n")
        if cfg == df.cfg.unique()[0]:
            L.append("\n### Paired Wilcoxon, krig_cortex vs spline (per output set and band)\n")
            L.append("| outputs | band | median dr | wins | p |\n|---|---|---|---|---|")
            for o in OUTPUT_SETS:
                for b in d.band.unique():
                    a = d[(d.outputs == o) & (d.band == b) & (d.method == "krig_cortex")].set_index("subject").r
                    s = d[(d.outputs == o) & (d.band == b) & (d.method == "spline")].set_index("subject").r
                    dl = (a - s).dropna()
                    p = wilcoxon(dl).pvalue if len(dl) > 5 and np.any(dl != 0) else np.nan
                    L.append(f"| {o} | {b} | {dl.median():+.3f} | {(dl > 0).sum()}/{len(dl)} | {p:.2g} |")
    out.with_suffix(".md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
