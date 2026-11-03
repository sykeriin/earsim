# earsim

Synthesise labelled EEG for any wearable electrode layout (earbuds, cEEGrid-style
around-ear arrays, headbands) from existing high-density scalp recordings, and
choose where the electrodes should go.

Full plan, rationale, alternatives and timeline: [`docs/PLAN.md`](docs/PLAN.md).

## Core idea in one line

With a linear inverse (MNE / eLORETA) the pipeline scalp → cortical sources → new
electrodes collapses to one matrix `T = R_tgt · G_tgt · K_scalp` (rank ≤ number of
scalp channels). Synthesis is therefore physics-informed spatial extrapolation.
It must beat spherical-spline interpolation to be worth anything, and it cannot
recover activity the scalp cap never saw.

## Status

See [`docs/RESULTS.md`](docs/RESULTS.md) for findings. Headline: physics-based synthesis helps a little
inside a scalp cap and does not help on real ear-EEG; zero-shot attention decoding from synthetic ear data fails.

## Modules (week 1 table, now extended)

| Module | State |
|---|---|
| `earsim/head.py` | fsaverage template, BEM, cached forward solutions |
| `earsim/layouts.py` | 10-05 positions, cEEGrid-like arc, crude in-ear proxies (Tier A: on the BEM scalp, no ear canal), candidate grid around each ear |
| `earsim/transfer.py` | explicit inverse kernel, transfer matrix `T`, reference operators, disk cache |
| `earsim/baselines.py` | spherical spline (tunable α), nearest neighbour |
| `earsim/data.py` | PhysioNet EEGMMIDB loader |
| `scripts/01_leave_channels_out.py` | E1: predict held-out near-ear cap channels; SNR/α chosen on dev subjects |
| sensor model, decoders, placement optimiser | not started (weeks 4–10) |

## Setup

```bash
pip install -e ".[dev]"
python -m pytest -q tests          # first run downloads fsaverage (~200 MB)
python scripts/01_leave_channels_out.py --subjects 1-109 --dev 1-20   # E1
python scripts/02_cache_covs.py && python scripts/02_kriging_e1.py      # E1 v2 (kriging)
python scripts/00_download_earsaad.py                                  # ~14 GB, Zenodo 16536441
python scripts/03_earsaad_cov.py && python scripts/03_e2.py            # E2
python scripts/05_e3.py --ref cross                                    # E3 (needs preprocessedData.zip unpacked in data/earsaad)
python scripts/06_placement_real.py                                    # E6 on real ear data
```

Data and the cache go to `../data` by default; override with `EARSIM_DATA=/path`.
