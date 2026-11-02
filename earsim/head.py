"""Template head model (fsaverage) and forward/lead-field computation.

Tier A head model: the stock fsaverage 3-layer BEM. Electrodes in or near the ear
are projected onto the smooth BEM scalp by MNE, so ear-canal geometry is NOT
modelled. Tier B (ear-augmented mesh) is a separate milestone; see docs/PLAN.md.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import mne
import numpy as np

DATA_DIR = Path(os.environ.get("EARSIM_DATA", Path(__file__).resolve().parents[2] / "data"))
SUBJECTS_DIR = DATA_DIR / "subjects"


@dataclass(frozen=True)
class TemplateHead:
    subjects_dir: Path
    trans: str
    src: mne.SourceSpaces
    bem: mne.bem.ConductorModel
    fiducials: dict  # name -> (3,) array, MRI coords in metres
    scalp_rr: np.ndarray  # (n, 3) BEM outer-skin vertices, MRI coords
    scalp_nn: np.ndarray


@lru_cache(maxsize=2)
def load_template(spacing: str = "ico-5") -> TemplateHead:
    """Load fsaverage, downloading it on first use."""
    SUBJECTS_DIR.mkdir(parents=True, exist_ok=True)
    fs_dir = Path(mne.datasets.fetch_fsaverage(subjects_dir=SUBJECTS_DIR, verbose=False))
    bem_dir = fs_dir / "bem"
    src = mne.read_source_spaces(bem_dir / f"fsaverage-{spacing}-src.fif", verbose=False)
    bem = mne.read_bem_solution(bem_dir / "fsaverage-5120-5120-5120-bem-sol.fif", verbose=False)
    fids, _ = mne.io.read_fiducials(bem_dir / "fsaverage-fiducials.fif")
    names = {1: "lpa", 2: "nasion", 3: "rpa"}
    fiducials = {names[int(f["ident"])]: np.asarray(f["r"]) for f in fids}
    skin = [s for s in bem["surfs"] if s["id"] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD][0]
    return TemplateHead(
        subjects_dir=SUBJECTS_DIR,
        trans="fsaverage",
        src=src,
        bem=bem,
        fiducials=fiducials,
        scalp_rr=skin["rr"].copy(),
        scalp_nn=skin["nn"].copy(),
    )


def make_info(ch_pos: dict[str, np.ndarray], sfreq: float = 160.0) -> mne.Info:
    """EEG Info with electrode positions given in fsaverage MRI coordinates (metres)."""
    head = load_template()
    montage = mne.channels.make_dig_montage(
        ch_pos=dict(ch_pos),
        nasion=head.fiducials["nasion"],
        lpa=head.fiducials["lpa"],
        rpa=head.fiducials["rpa"],
        coord_frame="mri",
    )
    info = mne.create_info(list(ch_pos), sfreq, "eeg")
    info.set_montage(montage)
    return info


def forward(info: mne.Info, n_jobs: int = 1) -> mne.Forward:
    """Free-orientation forward in surface orientation (needed by depth-weighted inverses).

    Cached on disk by electrode names + positions.
    """
    from . import cache

    pos = {ch["ch_name"]: ch["loc"][:3] for ch in info["chs"]}
    f = cache.path("fwd", cache.key(pos, mindist=5.0, src="ico-5"), "fwd.fif")
    if f.exists():
        fwd = mne.read_forward_solution(f, verbose=False)
    else:
        head = load_template()
        fwd = mne.make_forward_solution(
            info, trans=head.trans, src=head.src, bem=head.bem, eeg=True, meg=False,
            mindist=5.0, n_jobs=n_jobs, verbose=False,
        )
        mne.write_forward_solution(f, fwd, overwrite=True, verbose=False)
    return mne.convert_forward_solution(fwd, surf_ori=True, use_cps=True, verbose=False)


def fixed(fwd: mne.Forward) -> mne.Forward:
    return mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True,
                                        use_cps=True, copy=True, verbose=False)


def leadfield(info: mne.Info) -> tuple[np.ndarray, mne.Forward]:
    """Fixed-orientation (cortical normal) lead field G (n_channels x n_sources).

    Potentials are reference-free (MNE computes them against infinity); apply a
    reference explicitly afterwards.
    """
    fwd = fixed(forward(info))
    return fwd["sol"]["data"], fwd
