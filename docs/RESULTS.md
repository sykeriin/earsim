# Results so far (2 Nov 2026)

Short version: the physics-based transfer is a modest win over spherical-spline
interpolation **within a scalp cap**, and no win at all on **real ear-EEG**. The plan's
headline claim (decoders trained only on synthetic ear-EEG work on real ear-EEG) is **not
supported** by these experiments. Details, caveats and what each number can and cannot
say are below. Everything is reproducible from `scripts/` (commands in the README).

## E1 / E1 v2: leave-channels-out inside a scalp cap (EEGMMIDB, 84 test subjects)

- Plain MNE (identity source prior) does **not** beat a spline whose regulariser is tuned per band
  (`results/e1_eegbci_wide.md`). The first E1 write-up ("MNE wins alpha/beta") came from an
  under-tuned spline.
- **Kriging with fitted parcel variances** (`earsim/kriging.py`, cortex only, 68 variances fitted on
  20 dev subjects) beats the tuned spline when the whole lateral-temporal region (10 channels within
  7 cm of each ear) is held out: better in 5/6 bands (median dr delta +0.095, theta +0.019,
  alpha +0.007, gamma +0.031, broadband +0.045), no band significantly worse after Holm
  (beta -0.007, p=0.08). On the easier 6-channel hold-out it is mixed (beta significantly worse).
  Files: `results/e1v2.md`.
- Adding analytic ocular and muscle source groups made it **worse** (see `krig_full` in that
  file): alpha 0.873 -> 0.766 and beta 0.665 -> 0.522 on the ring split, with amplitude
  overshoot. Cause not isolated; suspects are the infinite-medium dipole model and guessed
  positions. Not part of the default model.

## E2: predict real ear channels from scalp (Ear-SAAD, leave-one-subject-out)

11 of 15 subjects usable (29 scalp channels, those with a bad source channel are dropped).
Outputs are reference-free differences. Ear electrode positions are template approximations;
the dataset ships none (`earsim/earsaad.py`). `results/e2.md`, sweep in `results/e2_sweep.md`.

- Kriging is **not better than spline**: e.g. cEEGrid local differences, broadband median r
  0.372 (kriging) vs 0.375 (spline) vs 0.463 (regression trained on other subjects' paired ear data).
  Kriging is slightly but significantly worse in beta/gamma local differences (p=0.002, 0.010,
  uncorrected, n=11).
- Cross-ear differences are predicted well by every method (r ~0.8 broadband, ~0.87 alpha): that is
  large-scale field structure, which a spline already captures.
- In-ear within-ear differences are essentially unpredictable (r ~0.01 for all physics/spline
  methods; 0.09 even for the paired-data regression).
- Not rescued by geometry: across 9 cEEGrid arc placements kriging-minus-spline stays within about
  +-0.03 except one geometry where the spline itself is unusually poor.
- Plan bar for go/no-go 2 ("beats spline and closes >=50% of the spline-to-oracle gap"): **not met**.

## E3: zero-shot attention decoding (Ear-SAAD, 15 subjects, device-realistic cross-ear reference)

Decoder trained on synthetic ear signals made from the subject's own scalp EEG, tested on the real
ear signals of a held-out trial. Window accuracy, mean over subjects (`results/e3_cross_oracle.md`):

| decoder trained on | 10 s | 30 s | 60 s |
|---|---|---|---|
| real ear (same subject, other trials; ceiling) | 0.579 | 0.622 | 0.679 |
| scalp (reference, not an ear) | 0.676 | 0.783 | 0.850 |
| synthetic ear, kriging | 0.511 | 0.512 | 0.510 |
| synthetic ear, spline | 0.502 | 0.497 | 0.500 |
| synthetic ear, linear map learned on OTHER subjects' real paired ear data | 0.514 | 0.501 | 0.519 |
| synthetic ear, linear map learned on the SAME subject's real paired ear data (other trials) | 0.537 | 0.552 | 0.564 |
| synthetic -> synthetic (control: train and test both synthetic) | 0.639 | 0.713 | 0.788 |

- Synthetic-trained decoders are at chance on real ear signals; kriging and spline are not
  distinguishable.
- The control shows the synthetic signals do carry attention information; it just does not
  transfer to the real ear signals.
- Even maps fitted on **real paired ear data** fail, so the limit is not specific to the physics
  model: for this task and these signals, a linear function of scalp EEG does not reproduce the
  attention-relevant part of the ear recording.
- With the scalp-mean reference (not device-realistic: the decoder input then contains the mean scalp
  signal), all 15 subjects (`results/e3_scalp_mean.md`): synthetic-kriging 0.517 / 0.535 / 0.538 at
  10 / 30 / 60 s (vs chance, p=0.21 at 60 s), synthetic-spline 0.505 / 0.512 / 0.530, maps learned on real
  paired data 0.524-0.576, real-ear decoder 0.579 / 0.633 / 0.671. A 3-subject smoke test had suggested
  0.63 for the kriging decoder; that was small-sample noise. Conclusion unchanged: synthetic-trained
  decoders do not approach real-ear ones under either reference.

## E6 (real data): where to put 1-3 electrodes in one ear (Ear-SAAD, 15 subjects)

`scripts/06_placement_real.py`, `results/e6_real.md`. Candidates are bipolar pairs inside one ear's
device. Greedy forward selection on the other 14 subjects, scored on the held-out subject with a
subject-specific decoder (leave-one-trial-out), 30 s windows. Chance = 0.5.

- cEEGrid, 3 pairs: greedy 0.587 (left) / 0.623 (right) vs 0.567 / 0.582 for random pairs from the same
  strip (+0.020, p=0.022; +0.041, p=0.0084; Wilcoxon, uncorrected). Weak but consistently above
  random in 11-13 of 15 subjects. The 19-electrode two-ear decoder reaches 0.62 at 30 s (E3), so a
  well-placed 3-pair single strip is not clearly worse than the full array here.
- In-ear mould (6 electrodes/ear): 0.51-0.55 for any k, essentially chance; greedy search helps only
  for a single pair (+0.03 left, +0.02 right, uncorrected) and not for k=2-3.
- Selected cEEGrid pairs are widely spaced (e.g. cEL5-cEL8, cER3-cER6), not neighbours.
- Pairs containing a dead channel are scored as unusable for that subject (contact quality counts).
- 12 comparisons were tested and none survives a Bonferroni correction (0.05/12): treat as
  suggestive. The reported ranking stability (Spearman 0.92-0.98 between leave-one-subject-out
  rankings) is inflated because those rankings share 13 of 15 subjects.
- Absolute decoding is weak (mean leave-one-trial-out correlation of the best pair is about 0.02), so
  differences between pairs are small compared to the subject-to-subject spread.

## What these results do and do not show

- They do not show that ear-EEG decoding is impossible from scalp data in general: only linear scalp
  maps (physics, spline, regression), one dataset, 15 subjects, one task, 20 Hz 1-9 Hz data.
- Template head, approximate ear geometry, no in-canal model (Tier A). A Tier-B ear-detailed head
  was not attempted; in-ear within-ear differences would need it, but the failing oracle maps
  suggest the problem is not only geometry.
- Hyper-parameters (spline alpha, MNE SNR, ridge strength) were chosen without test labels.
- n=11 (E2) and n=15 (E3) are small; p-values are uncorrected where noted.

## Where this leaves the project

The plan's pivot condition applies: drop the claim that physics-based synthesis from scalp data
replaces ear-EEG training data. What remains useful and supported: the transfer/kriging code as a
within-cap interpolation tool, the Ear-SAAD pipeline (robust covariances, decoding, placement search on
real ear data), and these negative results as a benchmark.
