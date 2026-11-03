"""E6 (real data): which 1-3 bipolar pairs inside one ear decode attention best?

Devices (one ear each): cEEGrid left/right (9/10 electrodes), in-ear mould left/right (6).
For each held-out subject the pairs are chosen by greedy search on the OTHER 14 subjects,
then scored on the held-out subject with a subject-specific decoder (leave-one-trial-out).
Baselines: random pairs from the same device (mean of 20 draws) and the single pair that
is best on average over training subjects. Output: results/e6_real.{csv,md}.
"""
import argparse
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from earsim import aad, earsaad, placement

warnings.filterwarnings("ignore")
DEVICES = {"cEEGrid_L": earsaad.CEEGRID_L, "cEEGrid_R": earsaad.CEEGRID_R,
           "in_ear_L": earsaad.IN_EAR_L, "in_ear_R": earsaad.IN_EAR_R}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default="1-15")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--win", type=int, default=30)
    ap.add_argument("--out", default="results/e6_real")
    args = ap.parse_args()
    a, _, b = args.subjects.partition("-")
    subs = list(range(int(a), int(b or a) + 1))
    data = {s: aad.load_prep(s) for s in subs}
    names = data[subs[0]][0]
    rng = np.random.default_rng(0)
    rows, rank_rows, t0 = [], [], time.time()

    for dev, group in DEVICES.items():
        pairs, cols = placement.pair_columns(names, group)
        stats = {s: placement.subject_stats(data[s][1], data[s][2], data[s][3], cols) for s in subs}
        score = placement.ScoreCache(stats)
        P = len(pairs)
        single = np.array([[score((j,))[s] for s in subs] for j in range(P)])  # (P, S)
        print(f"{dev}: {P} pairs, single-pair stats ({time.time() - t0:.0f}s)", flush=True)
        # ranking stability across held-out subjects (Spearman of leave-one-out mean rankings)
        loo = [np.delete(single, i, axis=1).mean(1) for i in range(len(subs))]
        rho = np.mean([spearmanr(loo[i], loo[j])[0] for i in range(len(subs)) for j in range(i)])
        for j in np.argsort(-single.mean(1))[:5]:
            rank_rows.append({"device": dev, "pair": "-".join(pairs[j]), "mean_loto_r": single[j].mean()})
        for hi, h in enumerate(subs):
            chosen = placement.greedy_select(score, P, args.k, exclude=h)
            top1 = [int(np.argmax(np.delete(single, hi, axis=1).mean(1)))]
            _, _, _, tr_h = data[h][0], data[h][1], data[h][2], data[h][3]
            sc = stats[h].scale
            for k in range(1, args.k + 1):
                ok, tot = placement.window_accuracy_for_subset(data[h][1], data[h][2], data[h][3], cols,
                                                               sc, tuple(chosen[:k]), args.win)
                rnd = []
                for _ in range(20):
                    r = tuple(rng.choice(P, k, replace=False))
                    o, t = placement.window_accuracy_for_subset(data[h][1], data[h][2], data[h][3], cols,
                                                                sc, r, args.win)
                    rnd.append(o / t)
                rows.append({"device": dev, "subject": h, "k": k, "greedy_acc": ok / tot,
                             "random_acc": float(np.mean(rnd)),
                             "greedy_pairs": ";".join("-".join(pairs[j]) for j in chosen[:k]),
                             "rank_stability_rho": float(rho)})
        print(f"{dev}: done ({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    df.to_csv(out.with_suffix(".csv"), index=False)
    L = [f"# E6 (real data): placement search for attention decoding, window {args.win}s\n",
         "Pairs chosen on the other subjects (leave-one-subject-out), scored on the held-out subject.\n",
         "Accuracy is the fraction of correctly classified windows; chance = 0.5.\n"]
    g = df.groupby(["device", "k"])[["greedy_acc", "random_acc"]].mean().round(3)
    L.append(g.to_markdown() + "\n")
    L.append("\n| device | k | greedy - random | wins | p |\n|---|---|---|---|---|")
    for (dev, k), d in df.groupby(["device", "k"]):
        diff = d.greedy_acc - d.random_acc
        p = wilcoxon(diff).pvalue if np.any(diff != 0) else np.nan
        L.append(f"| {dev} | {k} | {diff.mean():+.3f} | {(diff > 0).sum()}/{len(diff)} | {p:.2g} |")
    L.append("\nRanking stability (mean pairwise Spearman of leave-one-subject-out single-pair rankings):\n")
    L.append(df.groupby("device").rank_stability_rho.first().round(2).to_markdown() + "\n")
    L.append("\nTop single pairs by mean leave-one-trial-out correlation (all subjects):\n")
    L.append(pd.DataFrame(rank_rows).round(3).to_markdown(index=False) + "\n")
    out.with_suffix(".md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
