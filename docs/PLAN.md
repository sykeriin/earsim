> **Status banner (3 Nov 2026).** This is the original plan, kept for reference. Experiments since then did
> not support its headline claim (zero-shot decoders from synthetic ear-EEG). See `RESULTS.md` for findings and
> `NEXT_STEPS.md` for the revised plan. Revision 1 below (kriging) was implemented; its ocular/muscle source
> groups made predictions worse and are off by default.

# earsim: Form-factor simulator for wearable EEG (full project plan)

Date written: Mon 2 Nov 2026. Author context: solo, part-time. The plan assumes about 12–15 h/week.

---

## 0. Context: why this exists

Consumer EEG is moving into earbuds, headphones and headbands. Each device puts its electrodes somewhere new: the concha, the canal, behind the ear, along a headband line. Decoders trained on 10-20 scalp caps don't transfer to these sites. Labelled data from the new sites doesn't exist until someone builds hardware and runs a study, and that takes months per iteration.

**Goal:** take existing *labelled, high-density scalp* recordings and synthesise what any target electrode layout would have recorded, labels kept. Then use the synthetic data to (a) train decoders before the hardware exists and (b) choose where the 2–4 electrodes go.

**What this is NOT:** a way to create information the scalp recording lacks. See §2.1. That fact shapes the whole design and the validation.

---

## 1. Where things stand (honest status)

- Environment: Python 3.13, MNE 1.13.2, nibabel, sklearn, LightGBM installed. 4 CPU, 15 GB RAM, no GPU. That is fine for everything except large deep-learning sweeps.
- Downloaded and verified: fsaverage (BEM, ico-5 source space, fiducials) and PhysioNet EEGMMIDB subject 1 (64 ch, 160 Hz).
- Verified: MNE's `standard_1005` montage is already in fsaverage MRI coordinates (its LPA/RPA match the fsaverage fiducials exactly), so no coregistration step is needed.
- Verified: the fsaverage BEM scalp has ~112 vertices within 40 mm of each pre-auricular point. That is enough for an around-ear candidate grid.
- Scaffold written in `/home/user/earsim/earsim/` (not committed, not run end-to-end): `head.py`, `layouts.py`, `transfer.py`, `baselines.py`, `data.py`.
- **Known bug:** `transfer.inverse_kernel` calls `mne.make_eeg_average_ref_proj`, which doesn't exist at the top level in MNE 1.13. Fix: build a `RawArray`/`EvokedArray` on `fwd["info"]` and call `.set_eeg_reference("average", projection=True)`, then take its `info`.
- **Persistence risk:** this container is ephemeral and the repo has no remote. Work must be pushed to a GitHub repo you own, or it will be lost.

---

## 2. Technical framing: three things the idea doc gets wrong or understates

### 2.1 With a linear inverse, the whole pipeline is a single matrix
For MNE / eLORETA:
`ear = R_tgt · G_tgt · K_scalp · scalp = T · scalp`, with `T` of shape (n_target × n_scalp) and **rank ≤ n_scalp**.

Consequences:
- "Physics-based synthesis" here is **physics-informed spatial extrapolation**. The head model supplies the prior that a spline lacks. The fair baselines are therefore (a) spherical-spline interpolation and (b), where paired data exists, a learned linear regression scalp→ear, which is an *oracle ceiling* that the physics model cannot be expected to beat.
- Brain activity that the scalp cap does not see cannot be recovered. That includes local temporal-lobe activity and anything that is mostly EMG near the ear. The ear sits below and outside the convex hull of most caps, so this is extrapolation. That is the hard case.
- Synthesis is nearly free: one matrix product per layout. So we can generate thousands of layouts, and the placement optimiser can work on rows of `T` directly instead of re-running physics.

### 2.2 "Source estimate → project back to the standard channels → check the reconstruction" proves nothing
Re-projecting onto the *same* channels used for the inverse returns about the input as λ→0. It measures regularisation, not validity. **Replace it with leave-channels-out:** hold out the most inferior-lateral cap channels and predict them from the rest.
- EEGMMIDB channels: T9, T10, TP7, TP8, FT7, FT8.
- ANPHY-Sleep (83 ch): its lowest temporal channels.

### 2.3 Not every inverse method can be re-projected
dSPM and sLORETA output noise-normalised *statistics*, not current density. Pushing them through a lead field is physically meaningless. Allowed methods:
- MNE (L2) and eLORETA: linear, one kernel.
- Sparse inverses (MxNE, γ-MAP): non-linear, run per window, slow. Use them only as a comparison on short epochs.

### Additional corrections
- **Sleep source data must be high-density.** Classic PSG sets (Sleep-EDF, SHHS) have 2–6 EEG channels, and source estimation from those is meaningless. Use **ANPHY-Sleep** (83 ch, 29 healthy adults, OSF, Sci Data 2024).
- **The 70.1% JTEHM number is not a portable pass bar.** It was measured on a specific cohort (older people) and a specific dataset. The baseline "scalp-trained model applied to real ear-EEG" must be **re-run on the same test set** we use.
- **Motor imagery from the ear should sit near chance**, because motor cortex is far from the ear. So EEGMMIDB is for pipeline development and as a *negative control*, not a headline result. Ear-appropriate tasks are sleep, auditory attention, ASSR, auditory P300 and alpha/drowsiness.
- **Ethics approval in weeks 3–5 is unrealistic** at most institutions (typically 1–3 months), and it requires an institutional affiliation plus an amplifier. The own recording moves from "primary validation set" to an **optional phase after open data is exhausted**.

---

## 3. Data (verified where marked ✔)

| Role | Dataset | Content | Status |
|---|---|---|---|
| Dev / negative control | PhysioNet EEGMMIDB (Schalk 2004) | 64 ch, 109 subjects, motor imagery + eyes open/closed | ✔ downloaded |
| **Fidelity + ranking validation (primary)** | **Ear-SAAD**, Zenodo 16536441 (Geirnaert, Kappel, Kidmose, Sci Rep 2025) | Simultaneous scalp + 19-ch around-ear + in-ear; selective auditory attention; BIDS; ~12.6 GB | ✔ exists. Still to confirm: scalp channel count, ear electrode coordinates, stimulus envelopes |
| Sleep source data | ANPHY-Sleep (Wei et al., Sci Data 2024; OSF r26fh) | 83-ch HD PSG with scored stages and electrode positions | ✔ exists; licence to check |
| Sleep zero-shot test | EESM19 / EESM23 (Aarhus); EESM23 = OpenNeuro ds005178, CC0 | Ear-EEG sleep, scored; 10 subjects, ~1000 h | ✔ exists. Still to confirm: which channels are ear vs scalp |
| Auditory-attention source data | KU Leuven AAD (Das et al., 64 ch, 16 subjects); DTU AAD (Fuglsang et al., 64 ch, 18 subjects) | 64-ch BioSemi + audio envelopes | Known public sets; licence to check |
| Artifact library | EEGdenoiseNet (Zhang et al. 2021: EMG and EOG segments); TUH Artifact Corpus (has "chew" labels, needs registration) | Artifact segments | To download |

Week-2 emails go to the Aarhus group (Kidmose/Kappel: ear-electrode coordinates on a template, EESM channel maps) and the Oldenburg group (Debener/Bleichner: cEEGrid geometry, paired data). These double as the first cold emails.

---

## 4. System design

```
           labelled HD scalp data (Raw + annotations)
                        │
            ┌───────────▼───────────┐
            │ data.py loaders        │  standard montage in fsaverage MRI frame
            └───────────┬───────────┘
                        │
 head.py ──► lead fields G (fixed-orientation, ico-5, cached)
                        │
 transfer.py: K = inverse kernel(scalp montage, method, λ, depth)
              T = R_tgt · G_tgt · K          (cached by hash of positions+params)
                        │
 layouts.py: target layouts (cEEGrid-like, in-ear proxy, headband, candidate grid,
             device coordinates from CAD/photogrammetry when available)
                        │
 sensor.py:  reference scheme (in-ear / cross-ear / bipolar) → dry-electrode noise
             (1/f + white, impedance-scaled) → EMG/chew/motion from artifact library
             → band-limit / ADC; calibrated to real ear-EEG PSDs (unlabelled)
                        │
             synthetic labelled Raw/Epochs for the target layout
                ┌───────┴────────┐
 decode.py: EEGNet (braindecode), band-power + LightGBM, task-specific (AAD: linear stimulus reconstruction)
                │
 placement.py: greedy forward selection (baseline) and Concrete/Gumbel top-k selector
               over candidate rows of T (or candidate bipolar pairs), with device masks
                │
 validate.py: leave-channels-out, paired-data fidelity, zero-shot accuracy, layout ranking
```

**Key design decisions:**
- **Linear-operator core.** Everything flows through `T`, so it is cheap, cacheable and testable: linearity, reference invariance and rank checks.
- **Explicit reference operators.** Ear-EEG is dominated by its reference choice (electrodes 1–3 cm apart give tiny bipolar signals). The reference is a first-class matrix, not an afterthought.
- **Head-model domain randomisation** stands in for individual MRIs, which open datasets don't have. Jitter skull conductivity (ratio 1/25–1/80), head scale (±8%) and electrode position (±5 mm), and generate an ensemble of `T`s.
- **Sensor model calibrated on unlabelled target data.** We fit noise and artifact parameters so that synthetic PSD distributions match real ear-EEG. This uses no labels, so the zero-shot claim stays honest.
- **Packaging:** Python package plus `scripts/` per experiment; dataclass/YAML configs; results as CSV/JSON in `results/`; pytest unit tests.
- **Reuse instead of rewrite:** MNE (forward, inverse, BEM, spline via `mne.channels.interpolation._make_interpolation_matrix`), braindecode (EEGNet), MOABB where a dataset is supported, LightGBM.

---

## 5. Head-model tiers

- **Tier A (now):** stock fsaverage 3-layer BEM. MNE projects ear sites onto the smooth scalp. This is probably adequate for around-ear/cEEGrid sites and wrong for in-canal sites.
- **Tier B (only if Tier A fails in-ear fidelity on Ear-SAAD):** an ear-detailed head.
  - Options: (1) refine the BEM scalp around the pinna and canal; this is hard because BEM needs closed, nested, well-conditioned surfaces. (2) FEM (SimNIBS / DUNEuro) on a head with ear anatomy, e.g. MIDA (still to check whether it has canal detail), following Kappel, Makeig & Kidmose 2019.
  - Decision rule: do it only if the in-ear vs around-ear fidelity gap is large *and* attributable to geometry, not noise.
- **Ablation:** a 3-shell spherical model, to show whether a realistic head matters at all near the ear.

---

## 6. Experiments and pass bars

| # | Experiment | Metric | Pass bar | Go/no-go |
|---|---|---|---|---|
| E1 | Leave-channels-out on EEGMMIDB (109 subj) and ANPHY: predict T9/T10/TP7/TP8/FT7/FT8 | Per-band (δ, θ, α, β) correlation; relative RMSE after optimal gain; λ and depth sweep | MNE/eLORETA > spherical spline, paired Wilcoxon p<0.01 | **Go/no-go 1.** If physics doesn't beat spline here, the premise is weak. Report it and pivot to "spline + sensor model" |
| E2 | Ear-SAAD: scalp → real around-ear and in-ear channels | Same metrics; also the oracle linear-regression ceiling (cross-validated across subjects) | Beats spline; closes ≥50% of the spline→oracle gap | **Go/no-go 2** |
| E3 | Auditory attention zero-shot: train on synthetic ear from KUL/DTU 64-ch, test on real Ear-SAAD ear channels | AAD accuracy at fixed window lengths | > scalp-model-on-ear baseline; report gap to model trained on real ear | – |
| E4 | Sleep zero-shot: ANPHY → synthetic ear → test EESM | 5-class accuracy, κ | > re-run scalp→ear transfer baseline on the *same* EESM split | – |
| E5 | Alpha (eyes closed vs open), EEGMMIDB runs 1/2 → synthetic ear; test on any real ear set with rest blocks | AUC | sanity check | – |
| E6 | Placement: rank sub-layouts (subsets / bipolar pairs of Ear-SAAD's 19 around-ear channels) on synthetic vs real data | Spearman ρ between rankings | ρ > 0.6 | – |
| E7 | Ablations: no sensor model; no domain randomisation; spherical vs BEM; MNE vs eLORETA vs sparse | Δ on E2–E4 | descriptive | – |
| E8 | Negative control: motor imagery on ear layouts | Accuracy | ≈ chance (if not, suspect artifacts or leakage) | – |

Subject-wise splits everywhere. No hyperparameter tuning on test subjects. λ is chosen on E1 dev subjects only.

---

## 7. Placement optimiser

- Candidates: about 30–60 scalp sites per ear (`layouts.ear_candidates`) plus headband-line sites. Device masks enforce physical feasibility, e.g. one earbud = only the left-ear region and at most 3 contacts.
- Earbuds have no distant reference, so we select **bipolar pairs / contact sets**, not single electrodes.
- **Baseline:** greedy forward selection with a cheap decoder (band-power + LDA/LightGBM). For k ≤ 4 from ~60 candidates this is feasible and strong.
- **Research contribution:** a Concrete / Gumbel-softmax top-k selector over rows of `T_cand`, trained end-to-end with the decoder (Strypsteen & Bertrand 2021; Abid et al. 2019; Xie & Ermon 2019).
- Honest expectation: the differentiable selector may not beat greedy at k ≤ 4. Report the comparison either way.
- Output: the best layouts plus a robustness band across the head-model ensemble.

---

## 8. Timeline (start Mon 2 Nov 2026, ~12–15 h/week)

The idea doc says 10 weeks. With the late-December holidays and data wrangling, **12 weeks (to about 31 Jan 2027) is realistic**. Ten weeks is optimistic.

| Week | Dates | Work | Gate |
|---|---|---|---|
| 1 | Nov 2–8 | Fix scaffold bug; tests; E1 on EEGMMIDB (all subjects); λ/depth sweep | Go/no-go 1 |
| 2 | Nov 9–15 | Download Ear-SAAD, EESM23, ANPHY; build channel/coordinate maps; send emails | Data confirmed |
| 3 | Nov 16–22 | E2 fidelity on Ear-SAAD; E1 on ANPHY | Go/no-go 2 |
| 4 | Nov 23–29 | `sensor.py` v1; calibrate it against real ear PSDs; head-model domain randomisation | – |
| 5–6 | Nov 30–Dec 13 | `decode.py`; E3 auditory attention zero-shot; E5 alpha | – |
| 7–8 | Dec 14–27 (holiday slack) | E4 sleep zero-shot; re-run scalp-transfer baseline | – |
| 9–10 | Jan 4–17 | `placement.py`: greedy + Concrete; E6 ranking | – |
| 11 | Jan 18–24 | E7 ablations, E8 negative control; decide on Tier B | – |
| 12 | Jan 25–31 | Toolkit clean-up, docs, arXiv write-up, release | – |
| Optional | Feb 2027 → | Tier B head model; own paired recording (only with ethics approval and hardware) | – |

---

## 9. Alternatives considered and why they were not chosen

| Alternative | Why not (as the main approach) | Where it still appears |
|---|---|---|
| Deep scalp→ear generators (ACM 2025: CNN/LSTM/MAE) | Need paired data per subject and per layout, so they can't produce a layout that has never been recorded | Comparison on Ear-SAAD if code is available |
| Learned linear regression scalp→ear | Same limitation: it needs paired data for that exact layout | Oracle ceiling in E2 |
| Channel-agnostic foundation models (LaBraM, BIOT, etc.) | Still need some real data from the new layout to fine-tune, and they don't address placement | Future work: pre-train on synthetic layouts |
| Fully synthetic sources (SEREEGA, neural mass models) | No realistic task labels or between-subject variability | No |
| Individual-MRI head models | Open datasets don't have MRIs | Replaced by head-model domain randomisation |
| Spherical head | Fits the inferior and lateral head poorly, which is exactly where the ear is | Ablation |
| Sparse inverse as default | Non-linear, slow, unstable on continuous data; loses the single-matrix property | Comparison in E7 |
| dSPM / sLORETA | Statistics, not currents, so they can't be re-projected | Excluded |
| MATLAB (FieldTrip / Brainstorm) | Splits the stack away from PyTorch and braindecode | No |
| Own recording as primary validation | Ethics timeline, hardware, affiliation; Ear-SAAD already gives paired data | Optional phase |
| Evaluating on MI accuracy | The ear is far from motor cortex | Negative control |

---

## 10. Risks

1. **Extrapolation fails** (E1 or E2 no better than spline). Mitigation: it is cheap to find out in weeks 1–3, and a negative result is publishable as a benchmark. Pivot: spline + sensor model + placement.
2. **Real ear-EEG is dominated by noise/EMG**, so the sensor model decides the outcome. Mitigation: calibrate it on real ear PSDs; ablation E7.
3. **No ear electrode coordinates on the template.** Mitigation: emails; photo-based approximation; position-jitter randomisation.
4. **Dataset licences or access.** Check in week 2. EESM23 is CC0. Ear-SAAD and ANPHY terms still to read.
5. **Compute (no GPU).** EEGNet on ~100 subjects runs fine on CPU at 128–250 Hz. Heavy sweeps would need Colab or a rented GPU.
6. **Scope creep (Tier B mesh work).** It is gated behind the E2 results.

---

## 11. Related work to read (core)

- **Ear-EEG:**
  - Looney et al. 2012, *IEEE Pulse* (in-the-ear concept)
  - Kidmose et al. 2013, *IEEE TBME* (evoked potentials from ear-EEG)
  - Mikkelsen et al. 2019, *Sci Rep* (dry-contact ear-EEG sleep)
  - Kappel et al. 2019, *IEEE TBME* (dry-contact ear-EEG)
  - Kappel, Makeig & Kidmose 2019, *Front Neurosci* (ear-EEG forward models)
  - Geirnaert, Kappel & Kidmose 2025, *Sci Rep* (Ear-SAAD)
- **cEEGrid:** Debener et al. 2015, *Sci Rep*; Bleichner & Debener 2017, *Front Hum Neurosci*
- **Transfer:** JTEHM 2024 (scalp→ear sleep); ACM 2025 (ear-EEG synthesis); knowledge distillation scalp→ear sleep (arXiv 2211.02638)
- **Inverse / forward:**
  - Hämäläinen & Ilmoniemi 1994 (MNE)
  - Pascual-Marqui 2002/2007 (sLORETA / eLORETA)
  - Gramfort et al. 2012/2013 (MxNE; MNE-Python)
  - Wipf & Nagarajan 2009 (γ-MAP)
  - Vorwerk et al. 2014 (*NeuroImage*, head-model guidelines)
  - Huang, Parra & Haufe 2016 (New York Head)
- **Interpolation:** Perrin et al. 1989 (spherical splines)
- **Simulation:** Krol et al. 2018 (SEREEGA); Haufe & Ewald 2019 (simulation benchmark)
- **Decoders:** Lawhern et al. 2018 (EEGNet); Schirrmeister et al. 2017
- **AAD:** Das et al. 2016 (KUL dataset); Fuglsang et al. 2017 (DTU dataset)
- **Channel selection:** Strypsteen & Bertrand 2021 (*J Neural Eng*, Gumbel-softmax); Abid et al. 2019 (Concrete autoencoders); Xie & Ermon 2019
- **Artifacts:** Zhang et al. 2021 (EEGdenoiseNet); TUH artifact corpus
- **Data:** Schalk et al. 2004 (BCI2000 / EEGMMIDB); Wei et al. 2024 (ANPHY-Sleep)
- Not verifiable by me: the arXiv 2608.26998 virtual-iEEG survey is past my knowledge cutoff, so read it yourself before citing it.

---

## 12. Verification of the build (how we know each piece works)

- Unit tests (pytest):
  - `T` is linear.
  - Adding a constant to all scalp channels leaves the output unchanged (average reference).
  - rank(T) ≤ n_scalp.
  - Bipolar operator rows sum to 0.
  - Re-projecting onto the source montage with small λ ≈ identity on avg-ref data. This is a sanity check only.
  - Layouts lie on the scalp surface (distance < 1 mm).
- Week-1 end-to-end: run `python scripts/01_leave_channels_out.py --subjects 1-109`. It writes `results/e1_eegbci.csv` plus a per-band summary table of MNE vs eLORETA vs spline vs nearest-neighbour, with paired statistics.
- Every experiment script writes its config plus git hash next to its results.

## 13. First build steps after approval (week 1)

1. Fix the average-reference projection bug in `earsim/transfer.py`.
2. Add a `cache.py` keyed on positions + params.
3. Add `tests/test_transfer.py` covering the checks in §12.
4. Add `scripts/01_leave_channels_out.py` (EEGMMIDB, λ/depth sweep, per-band metrics, Wilcoxon).
5. Write `docs/PLAN.md` (this plan) and a README into the repo.
6. Commit; push once you provide a GitHub repo.
