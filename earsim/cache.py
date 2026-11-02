"""Disk cache for expensive, deterministic artefacts (forwards, transfer matrices).

Keys are hashes of electrode names + positions (rounded to 0.1 mm) + parameters,
so identical montages across subjects reuse the same matrices.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .head import DATA_DIR

CACHE_DIR = DATA_DIR / "cache"


def key(pos: dict[str, np.ndarray] | None = None, **params) -> str:
    h = hashlib.sha1()
    if pos is not None:
        for name, p in pos.items():  # order matters: rows of the matrices follow it
            h.update(name.encode())
            h.update(np.round(np.asarray(p, float) * 1e4).astype(np.int64).tobytes())
    h.update(json.dumps(params, sort_keys=True, default=str).encode())
    return h.hexdigest()[:16]


def path(kind: str, k: str, ext: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{kind}-{k}.{ext}"
