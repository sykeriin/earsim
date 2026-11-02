"""Robust per-subject band covariances for Ear-SAAD (scalp + in-ear + around-ear).

The raw recordings contain dead channels (flat), channels with huge artefact power
(e.g. a cEEGrid electrode at 3 mV) and many artefact windows. For each subject:
  1. resample to 250 Hz, drop EOG,
  2. mark bad channels per modality group (flat, or robust z of log median power > 3.5),
  3. mark bad 2 s windows (>= 10 % of good channels with robust z > 4 against their own
     window-power distribution) on the 1-45 Hz signal,
  4. average x x^T over good windows, per band.
Output data/cache/cov-earsaad-sub{NN}.npz: C[band, 64, 64] (NaN for bad channels),
good[64], n_windows_kept/total, ch_names, bands. Subjects whose download may still be in
progress are skipped unless --force.
"""
import argparse
import time
from pathlib import Path

import mne
import numpy as np
from scipy.signal import butter, sosfiltfilt

from earsim import cache
from earsim.earsaad import GROUPS
from earsim.head import DATA_DIR

mne.set_log_level("ERROR")
BANDS = {"delta": (1, 4), "theta": (4, 8), "alpha": (8, 13), "beta": (13, 30),
         "gamma": (30, 45), "broad": (1, 45)}
SF = 250
WIN = 2 * SF
ROOT = DATA_DIR / "earsaad/bids_dataset"


def robust_z(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826 + 1e-12
    return (x - med) / mad


def process(sub: int, force=False):
    d = ROOT / f"sub-{sub:02d}/ses-01/eeg"
    f = d / f"sub-{sub:02d}_ses-01_task-selectiveAttention_eeg.set"
    raw = mne.io.read_raw_eeglab(f, preload=True)
    names = raw.ch_names
    raw.resample(SF)
    X = raw.get_data() * 1e6  # uV
    X -= X.mean(1, keepdims=True)
    n_ch, n_t = X.shape
    nw = n_t // WIN
    # broad-band windows for QC
    sos = butter(4, [1, 45], "bandpass", fs=SF, output="sos")
    Xb = sosfiltfilt(sos, X)
    Wv = Xb[:, : nw * WIN].reshape(n_ch, nw, WIN).var(2) + 1e-12  # (ch, win)
    good = np.zeros(n_ch, bool)
    for grp, chs in GROUPS.items():
        idx = [names.index(c) for c in chs]
        med = np.median(Wv[idx], axis=1)
        flat = med < 1e-6
        # bad = flat, or power > 16x (amplitude > 4x) the group median. Quiet channels are
        # kept: those next to the common reference (Fpz) are legitimately low-amplitude.
        ok = ~flat & (med < 16 * np.median(med[~flat]))
        good[idx] = ok
    lw = np.log(Wv)
    zw = np.stack([robust_z(lw[c]) for c in range(n_ch)])
    frac = (zw[good] > 4).mean(0)
    keep = frac < 0.10
    cov = np.full((len(BANDS), n_ch, n_ch), np.nan)
    for bi, (lo, hi) in enumerate(BANDS.values()):
        sos = butter(4, [lo, min(hi, SF / 2 - 1)], "bandpass", fs=SF, output="sos")
        Xf = sosfiltfilt(sos, X[good])
        Xw = Xf[:, : nw * WIN].reshape(good.sum(), nw, WIN)[:, keep].reshape(good.sum(), -1)
        c = Xw @ Xw.T / Xw.shape[1]
        gi = np.where(good)[0]
        cov[bi][np.ix_(gi, gi)] = c
    p = cache.path("cov", f"earsaad-sub{sub:02d}", "npz")
    np.savez(p, C=cov, good=good, kept=int(keep.sum()), total=int(nw), ch_names=np.array(names),
             bands=np.array(list(BANDS)))
    bad = [n for n, g in zip(names, good) if not g]
    print(f"sub {sub:02d}: windows kept {keep.sum()}/{nw}; bad channels {bad}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default="1-15")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    a, _, b = args.subjects.partition("-")
    subs = list(range(int(a), int(b or a) + 1))
    t0 = time.time()
    present = [s for s in subs if (ROOT / f"sub-{s:02d}/ses-01/eeg").exists()]
    for s in present:
        out = cache.path("cov", f"earsaad-sub{s:02d}", "npz")
        if out.exists() and not args.force:
            continue
        later = any(t > s for t in present)
        if not later and "done" not in Path("results/earsaad_download.log").read_text():
            print(f"sub {s:02d}: download may be in progress, skipping")
            continue
        process(s)
        print(f"  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
