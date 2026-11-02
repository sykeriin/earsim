"""Dataset loaders. Each returns an mne.io.Raw with a montage in fsaverage MRI coords."""
from __future__ import annotations

import warnings

import mne

from .head import DATA_DIR


def load_eegbci(subject: int, runs=(4, 8, 12)) -> mne.io.BaseRaw:
    """PhysioNet EEG Motor Movement/Imagery (Schalk 2004), 64 ch, 160 Hz.

    Runs 4/8/12: imagined left vs right fist. Run 1/2: eyes open/closed baseline.
    """
    from mne.datasets import eegbci
    files = eegbci.load_data(subject, list(runs), path=DATA_DIR, update_path=False, verbose=False)
    raws = [mne.io.read_raw_edf(f, preload=True, verbose=False) for f in files]
    raw = mne.concatenate_raws(raws, verbose=False)
    eegbci.standardize(raw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        raw.set_montage("standard_1005", verbose=False)
    raw.annotations.rename({"T0": "rest", "T1": "left", "T2": "right"})
    return raw
