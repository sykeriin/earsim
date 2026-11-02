"""E1: leave-channels-out validation on PhysioNet EEGMMIDB.

Hold out the most inferior-lateral cap channels (closest to the ear) and predict
them from the remaining channels with
  - source transfer (MNE, eLORETA) through the fsaverage BEM,
  - spherical-spline interpolation (Perrin 1989),
  - nearest neighbour.
Ground truth is the recorded held-out signal, both referenced to the mean of the
kept channels. SNR (lambda) is selected on dev subjects only; stats on the rest.

Usage:
  python scripts/01_leave_channels_out.py --subjects 1-109 --dev 1-20
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import mne
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from scipy.stats import wilcoxon

from earsim.baselines import nearest_neighbour, spherical_spline
from earsim.data import load_eegbci
from earsim.transfer import build_transfer

mne.set_log_level("ERROR")

BANDS = {"delta": (1, 4), "theta": (4, 8), "alpha": (8, 13), "beta": (13, 30),
         "gamma": (30, 45), "broad": (1, 45)}
DEFAULT_HELD = ["T9", "T10", "TP7", "TP8", "FT7", "FT8"]
# Commonly excluded in MOABB/braindecode: different sfreq or broken annotations.
BAD_SUBJECTS = {88, 89, 92, 100, 104}
METHODS = [("MNE", 1.0), ("MNE", 3.0), ("MNE", 10.0), ("eLORETA", 3.0)]
SPLINE_ALPHAS = [1e-5, 1e-7, 1e-9]  # 1e-5 is MNE's interpolate_bads default


def parse_range(s: str) -> list[int]:
    out = []
    for part in s.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def metrics(pred: np.ndarray, true: np.ndarray) -> dict[str, np.ndarray]:
    pc = pred - pred.mean(1, keepdims=True)
    tc = true - true.mean(1, keepdims=True)
    r = (pc * tc).sum(1) / np.sqrt((pc**2).sum(1) * (tc**2).sum(1))
    rel_rmse = np.sqrt(((pc - tc) ** 2).mean(1)) / tc.std(1)
    gain = pc.std(1) / tc.std(1)  # <1 means the method under-predicts amplitude
    return {"r": r, "rel_rmse": rel_rmse, "gain": gain}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default="1-109")
    ap.add_argument("--dev", default="1-20", help="subjects used only for SNR selection")
    ap.add_argument("--runs", default="1,2,4")
    ap.add_argument("--held", default=",".join(DEFAULT_HELD))
    ap.add_argument("--out", default="results/e1_eegbci")
    args = ap.parse_args()

    subjects = [s for s in parse_range(args.subjects) if s not in BAD_SUBJECTS]
    dev = set(parse_range(args.dev))
    runs = tuple(int(r) for r in args.runs.split(","))
    held = args.held.split(",")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    raw0 = load_eegbci(subjects[0], runs)
    pos = raw0.get_montage().get_positions()["ch_pos"]
    kept = {k: v for k, v in pos.items() if k not in held}
    tgt = {k: pos[k] for k in held}

    t0 = time.time()
    ops = {}
    for m, snr in METHODS:
        ops[f"{m}_snr{snr:g}"] = build_transfer(kept, tgt, method=m, snr=snr).T
        print(f"built {m} snr={snr} ({time.time() - t0:.0f}s)", flush=True)
    for a in SPLINE_ALPHAS:
        ops[f"spline_a{a:g}"] = spherical_spline(kept, tgt, alpha=a)
    ops["nearest"] = nearest_neighbour(kept, tgt)

    rows = []
    for s in subjects:
        try:
            raw = raw0 if s == subjects[0] else load_eegbci(s, runs)
        except Exception as e:  # noqa: BLE001
            print(f"subject {s}: skipped ({e})")
            continue
        sf = raw.info["sfreq"]
        X = raw.get_data(picks=list(kept))
        Y = raw.get_data(picks=held)
        ref = X.mean(0, keepdims=True)
        X, Y = X - ref, Y - ref
        for band, (lo, hi) in BANDS.items():
            sos = butter(4, [lo, min(hi, sf / 2 - 1)], "bandpass", fs=sf, output="sos")
            Xb, Yb = sosfiltfilt(sos, X), sosfiltfilt(sos, Y)
            for name, W in ops.items():
                met = metrics(W @ Xb, Yb)
                for i, ch in enumerate(held):
                    rows.append({"subject": s, "split": "dev" if s in dev else "test",
                                 "band": band, "method": name, "channel": ch,
                                 **{k: float(v[i]) for k, v in met.items()}})
        print(f"subject {s} done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out.with_suffix(".csv"), index=False)
    summarise(df, out)
    git = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                         text=True).stdout.strip() or "uncommitted"
    json.dump({"args": vars(args), "methods": METHODS, "held": held, "git": git,
               "n_subjects": int(df.subject.nunique())},
              open(out.with_suffix(".json"), "w"), indent=2)


def summarise(df: pd.DataFrame, out: Path) -> None:
    subj = df.groupby(["split", "subject", "band", "method"])[["r", "rel_rmse", "gain"]] \
             .mean().reset_index()
    dev = subj[subj.split == "dev"]
    mne_methods = [m for m in subj.method.unique() if m.startswith("MNE_")]
    best = dev[dev.method.isin(mne_methods)].groupby("method").r.mean().idxmax() \
        if len(dev) else "MNE_snr3"
    spl = [m for m in subj.method.unique() if m.startswith("spline_")]
    best_spl = dev[dev.method.isin(spl)].groupby("method").r.mean().idxmax() \
        if len(dev) else "spline_a1e-05"
    test = subj[subj.split == "test"] if (subj.split == "test").any() else subj
    lines = [f"# E1 leave-channels-out (EEGMMIDB)\n",
             f"Subjects: {test.subject.nunique()} test, {dev.subject.nunique()} dev. "
             f"Selected on dev: MNE **{best}**, spline **{best_spl}**.\n",
             "Median over test subjects of the channel-mean correlation (r), relative RMSE "
             "and gain (pred std / true std).\n"]
    tab = test.groupby(["band", "method"])[["r", "rel_rmse", "gain"]].median().round(3)
    lines.append(tab.unstack("method")["r"].to_markdown())
    lines.append("\n\nRelative RMSE:\n")
    lines.append(tab.unstack("method")["rel_rmse"].to_markdown())
    lines.append("\n\nGain:\n")
    lines.append(tab.unstack("method")["gain"].to_markdown())
    lines.append("\n\n## Paired Wilcoxon on r (test subjects)\n")
    lines.append("| band | comparison | median Δr | wins | p |\n|---|---|---|---|---|")
    for band in test.band.unique():
        piv = test[test.band == band].pivot(index="subject", columns="method", values="r")
        for a in (best, "eLORETA_snr3"):
            if a not in piv or best_spl not in piv:
                continue
            d = (piv[a] - piv[best_spl]).dropna()
            p = wilcoxon(d).pvalue if len(d) > 5 and np.any(d != 0) else np.nan
            lines.append(f"| {band} | {a} − {best_spl} | {d.median():+.4f} | "
                         f"{(d > 0).sum()}/{len(d)} | {p:.2g} |")
    out.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
