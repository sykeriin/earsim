# Next steps

Status as of 3 Nov 2026: see [`RESULTS.md`](RESULTS.md). Short version: physics-based synthesis helps a
little inside a scalp cap, does not beat spline interpolation on real ear channels, and zero-shot attention
decoding from synthetic ear data is at chance (even with maps learned from real paired data). The original
plan ([`PLAN.md`](PLAN.md)) is kept for reference; its headline claim is not supported yet.

## 1. Decisive experiment: sleep staging (E4), the most favourable remaining test

Attention decoding depends on small, fast, local signals that scalp data do not predict. Sleep staging
depends on large, slow, spatially smooth activity, which is exactly what scalp-to-ear transfer predicts well
(delta and cross-ear differences reached r about 0.8 in E2). If synthetic ear data cannot train a sleep
stager that works on real ear-EEG, the idea is in serious trouble; if it can, the project has a defensible
niche.

- Data: ANPHY-Sleep (83-channel scalp PSG, 29 subjects, OSF r26fh; licence to check) as the source;
  EESM19/EESM23 (OpenNeuro ds005178, CC0; first confirm which channels are ear and which scalp) as the
  real ear test. Re-run the scalp-trained-model-on-real-ear baseline on the same split (do not reuse the
  70.1 % figure from the literature).
- Reuse: `earsim/kriging.py` (fit parcel variances on ANPHY covariances, then predict ear sites), the
  ear-position helpers in `earsim/earsaad.py` (add EESM electrode sites), `earsim/aad.py`-style
  leave-one-subject-out scaffolding.
- Pass bar: synthetic-trained stager beats the scalp-model-on-ear baseline and a spline-synthetic
  baseline, paired across subjects.
- Decision rule: fail -> stop developing synthesis and publish the benchmark (section 4); pass ->
  continue to the sensor/artifact model and placement optimiser (original plan, weeks 4-10).
- Effort: about 1-2 weeks, mostly data wrangling.

## 2. Close known gaps in the existing experiments (1-3 days)

- E2 uses 11 of 15 subjects because every source channel must be good; support per-subject source sets.
- The covariance-based E1 metric reproduces the earlier time-series MNE number to 0.003, not exactly; find
  out why (cached transfer order or time-mean handling are the candidates).
- Add the planned test that kriging with an identity source prior reproduces MNE.
- Apply a multiple-comparison correction to the E6 placement tests (currently uncorrected; none passes
  Bonferroni).
- Try non-linear decoders (EEGNet) on E3; a linear map may fail where a network on the same inputs does not.
  Expect low odds given the paired-data oracle result, but it is cheap to rule out.
- Re-run E3 `--ref local` on all subjects (only cross and scalp_mean are in the repo).
- Re-check any conclusion drawn from fewer than about 10 subjects: the 3-subject scalp-mean smoke test
  was misleading (see RESULTS.md).

## 3. Real ear geometry

The dataset ships no electrode coordinates, so every ear-site result rests on guessed positions.
- Email the Aarhus group (Kidmose, Kappel) for template coordinates of the in-ear and cEEGrid sites.
- Only after coordinates exist, and only if the in-ear gap is attributable to geometry, consider an
  ear-detailed head (FEM) as in the original plan's Tier B. The failing paired-data oracle suggests geometry
  is not the main bottleneck, so this is lower priority than section 1.

## 4. Benchmark / negative-result write-up

If section 1 does not rescue the approach, the supportable contribution is: an open benchmark and toolkit
showing that scalp-to-ear synthesis (physics, spline, learned linear) does not transfer to attention decoding
on real ear-EEG, plus a real-data electrode-placement analysis. Outline: problem, toolkit, E1-E3 and E6
results, limitations, reproduction commands. Draft in `docs/RESULTS.md`; add figures (the repo has none yet).

## 5. Not recommended for now

- More physics (ocular/muscle sources, ear canal, domain randomisation): no evidence yet that physics is the
  bottleneck, and the first attempt at extra sources made predictions worse.
- Recording a new paired dataset: needs ethics approval, an affiliation and an amplifier; do it only if
  section 1 is positive and public data runs out.
- Differentiable placement optimiser (Concrete/Gumbel): greedy search is already cheap and the effects being
  optimised are small; revisit only if synthetic data become informative.

## Reproducing and persisting

Raw data and caches are not in the repo (about 15 GB). Commands are in the README; the container used so far
is ephemeral, so keep `data/` on your own disk (set `EARSIM_DATA`) or re-download with
`scripts/00_download_earsaad.py` plus the `preprocessedData.zip` from Zenodo record 16536441.
