"""Cache per-subject, per-band 64x64 channel covariances for EEGMMIDB.

Every E1-style metric and the kriging fit need only second-order statistics, so one
pass over the raw files is enough. Output: data/cache/cov_eegbci.npz with
  C[subject_index, band, 64, 64], subjects, bands, ch_names
"""
import argparse
import time

import mne
import numpy as np
from scipy.signal import butter, sosfiltfilt

from earsim import cache
from earsim.data import load_eegbci

mne.set_log_level("ERROR")
BANDS = {"delta": (1, 4), "theta": (4, 8), "alpha": (8, 13), "beta": (13, 30),
         "gamma": (30, 45), "broad": (1, 45)}
BAD = {88, 89, 92, 100, 104}

ap = argparse.ArgumentParser()
ap.add_argument("--subjects", default="1-109")
ap.add_argument("--runs", default="1,2,4")
args = ap.parse_args()
a, _, b = args.subjects.partition("-")
subjects = [s for s in range(int(a), int(b or a) + 1) if s not in BAD]
runs = tuple(int(r) for r in args.runs.split(","))

out, kept, names, t0 = [], [], None, time.time()
for s in subjects:
    try:
        raw = load_eegbci(s, runs)
    except Exception as e:  # noqa: BLE001
        print(f"subject {s}: skipped ({e})", flush=True)
        continue
    if names is None:
        names = raw.ch_names
    assert raw.ch_names == names
    sf = raw.info["sfreq"]
    X = raw.get_data()
    X -= X.mean(0, keepdims=True)
    Cs = []
    for lo, hi in BANDS.values():
        sos = butter(4, [lo, min(hi, sf / 2 - 1)], "bandpass", fs=sf, output="sos")
        Xb = sosfiltfilt(sos, X)
        Cs.append(Xb @ Xb.T / Xb.shape[1])
    out.append(np.stack(Cs))
    kept.append(s)
    print(f"subject {s} ({time.time() - t0:.0f}s)", flush=True)
f = cache.path("cov", "eegbci", "npz")
np.savez(f, C=np.stack(out), subjects=np.array(kept), bands=np.array(list(BANDS)),
         ch_names=np.array(names))
print("saved", f, np.stack(out).shape)
