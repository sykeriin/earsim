"""E1 v2: kriging with an augmented source space vs MNE vs tuned spline.

Two splits on EEGMMIDB (metrics computed from cached 64x64 band covariances):
  interp  hold out T9 T10 TP7 TP8 FT7 FT8 (as E1)
  ring    hold out every channel within 7 cm of either pre-auricular point (10 channels:
          the whole lateral-temporal region; real extrapolation)

For every method and band, hyper-parameters (MNE SNR, spline alpha) are chosen on dev
subjects (1-20) *per band*; the kriging variance parameters theta are fitted per band
on dev subjects' covariances of the kept channels only. Stats on test subjects only.

Usage: python scripts/02_kriging_e1.py [--splits interp,ring] [--out results/e1v2]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from earsim import kriging as kr
from earsim import layouts, sources
from earsim.baselines import nearest_neighbour, spherical_spline
from earsim.head import load_template
from earsim.transfer import build_transfer

warnings.filterwarnings("ignore")
INTERP = ["T9", "T10", "TP7", "TP8", "FT7", "FT8"]
SPLINE_ALPHAS = [1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2]
MNE_SNRS = [1.0, 3.0, 10.0]


def ring_channels(names, pos, radius=0.07):
    h = load_template()
    P = np.array([pos[n] for n in names])
    d = np.min([np.linalg.norm(P - h.fiducials[k], axis=1) for k in ("lpa", "rpa")], axis=0)
    return [n for n, x in zip(names, d) if x < radius]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", default="interp,ring")
    ap.add_argument("--dev-max", type=int, default=20)
    ap.add_argument("--out", default="results/e1v2")
    ap.add_argument("--ridge", type=float, default=1e-3)
    ap.add_argument("--no-extra", action="store_true", help="skip the ocular/muscle ablation row")
    args = ap.parse_args()

    z = np.load(Path(__import__("earsim.head", fromlist=["DATA_DIR"]).DATA_DIR) / "cache/cov-eegbci.npz")
    C, subjects, bands, names = z["C"], z["subjects"], list(z["bands"]), [str(n) for n in z["ch_names"]]
    pos = layouts.standard(names)
    dev_mask = subjects <= args.dev_max
    print(f"{len(subjects)} subjects ({dev_mask.sum()} dev), bands {bands}", flush=True)

    t0 = time.time()
    n_ctx = None
    bnames, S_all = sources.all_bases(pos, use_extra=True)
    n_ctx = sum(n.startswith("ctx:") for n in bnames)
    print(f"bases: {n_ctx} cortical + {len(bnames) - n_ctx} extra ({time.time() - t0:.0f}s)", flush=True)

    rows, fits = [], {}
    for split in args.splits.split(","):
        held = INTERP if split == "interp" else ring_channels(names, pos)
        kept = [n for n in names if n not in held]
        ki = np.array([names.index(n) for n in kept])
        ti = np.array([names.index(n) for n in held])
        print(f"\n== split {split}: hold out {held} ({len(kept)} kept)", flush=True)
        kpos = {n: pos[n] for n in kept}
        tpos = {n: pos[n] for n in held}

        ops = {}
        for snr in MNE_SNRS:
            ops[f"MNE|snr={snr:g}"] = build_transfer(kpos, tpos, method="MNE", snr=snr).T
        for a in SPLINE_ALPHAS:
            ops[f"spline|a={a:g}"] = spherical_spline(kpos, tpos, alpha=a)
        ops["nearest|-"] = nearest_neighbour(kpos, tpos)
        print(f"  baselines built ({time.time() - t0:.0f}s)", flush=True)

        B = kr.mean_zero_basis(len(kept))
        variants = {"krig_full": np.arange(len(bnames))}
        if not args.no_extra:
            variants["krig_cortex"] = np.arange(n_ctx)
        for bi, band in enumerate(bands):
            Ckk = np.stack([B.T @ C[s, bi][np.ix_(ki, ki)] @ B for s in range(len(subjects))])
            for vname, sel in variants.items():
                S_B = np.stack([B.T @ S_all[k][np.ix_(ki, ki)] @ B for k in sel])
                th = kr.fit_theta(S_B, Ckk[dev_mask], ridge=args.ridge)
                Sig = np.tensordot(th[:-1], S_all[sel], 1) + th[-1] * np.eye(len(names))
                ops[f"{vname}|{band}"] = kr.predict_operator(Sig, ki, ti)
                fits[f"{split}|{band}|{vname}"] = {"theta": th.tolist(),
                                                   "names": [bnames[k] for k in sel] + ["noise"]}
            print(f"  fitted {band} ({time.time() - t0:.0f}s)", flush=True)

        for si, s in enumerate(subjects):
            for bi, band in enumerate(bands):
                for key, W in ops.items():
                    fam, par = key.split("|")
                    if fam.startswith("krig") and par != band:
                        continue
                    m = kr.metrics_from_cov(W, C[si, bi], ki, ti)
                    rows.append({"split": split, "subject": int(s), "dev": bool(s <= args.dev_max),
                                 "band": band, "family": fam, "param": par if not fam.startswith("krig") else "-",
                                 **{k: float(np.mean(v)) for k, v in m.items()}})
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out.with_suffix(".csv"), index=False)
    json.dump(fits, open(out.with_suffix(".fits.json"), "w"))
    summarise(df, out)
    git = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump({"args": vars(args), "git": git, "ring": ring_channels(names, pos)},
              open(out.with_suffix(".json"), "w"), indent=2)


def summarise(df: pd.DataFrame, out: Path) -> None:
    L = ["# E1 v2: kriging vs MNE vs spline (EEGMMIDB, covariance-based metrics)\n",
         "Per band, MNE SNR and spline alpha are chosen on dev subjects (1-20). Statistics on test subjects.\n"]
    for split in df.split.unique():
        d = df[df.split == split]
        bands = list(dict.fromkeys(d.band))
        # pick best parameter per family and band on dev
        best = {}
        for fam in ("MNE", "spline"):
            for b in bands:
                dd = d[(d.family == fam) & (d.band == b) & d.dev]
                best[(fam, b)] = dd.groupby("param").r.mean().idxmax()
        def pick(fam, b):
            if fam in ("MNE", "spline"):
                return d[(d.family == fam) & (d.band == b) & (d.param == best[(fam, b)])]
            return d[(d.family == fam) & (d.band == b)]
        fams = ["nearest", "MNE", "spline"] + [f for f in d.family.unique() if f.startswith("krig")]
        L.append(f"\n## split: {split}\n")
        for metric in ("r", "rel_rmse", "gain"):
            tab = {}
            for fam in fams:
                tab[fam] = [pick(fam, b)[~pick(fam, b).dev][metric].median() for b in bands]
            L.append(f"\nMedian {metric} (test subjects)\n")
            L.append(pd.DataFrame(tab, index=bands).round(3).to_markdown() + "\n")
        L.append("\nChosen on dev: " + ", ".join(f"{fam}/{b}={best[(fam, b)]}"
                                                  for fam in ("MNE", "spline") for b in bands) + "\n")
        L.append("\n### Paired Wilcoxon on r vs tuned spline (test subjects; Holm over bands per method)\n")
        L.append("| method | band | median Δr | wins | p | Holm p |\n|---|---|---|---|---|---|")
        for fam in fams:
            if fam == "spline":
                continue
            ps, rec = [], []
            for b in bands:
                a, s_ = pick(fam, b), pick("spline", b)
                a = a[~a.dev].set_index("subject").r
                s_ = s_[~s_.dev].set_index("subject").r
                dlt = (a - s_).dropna()
                p = wilcoxon(dlt).pvalue if np.any(dlt != 0) else np.nan
                ps.append(p)
                rec.append((b, dlt.median(), int((dlt > 0).sum()), len(dlt), p))
            order = np.argsort(ps)
            holm = np.empty(len(ps))
            run = 0.0
            for rank, i in enumerate(order):
                run = max(run, min(1.0, (len(ps) - rank) * ps[i]))
                holm[i] = run
            for (b, md, w, n, p), h in zip(rec, holm):
                L.append(f"| {fam} | {b} | {md:+.4f} | {w}/{n} | {p:.2g} | {h:.2g} |")
    out.with_suffix(".md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
