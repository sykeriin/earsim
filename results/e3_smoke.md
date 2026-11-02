# E3: zero-shot attention decoding from synthetic ear-EEG (Ear-SAAD, 15 subjects)

Accuracy = fraction of windows where the attended stream is picked; chance = 0.5.

Per-subject accuracy, mean over subjects (and number of subjects above chance by binomial test p<0.05).

| method              |   win_s |   mean_acc |   median_acc | subjects>chance   |
|:--------------------|--------:|-----------:|-------------:|:------------------|
| real_ear            |      10 |      0.592 |        0.581 | 3/3               |
| real_ear            |      30 |      0.628 |        0.65  | 2/3               |
| real_ear            |      60 |      0.644 |        0.667 | 2/3               |
| scalp               |      10 |      0.715 |        0.717 | 3/3               |
| scalp               |      30 |      0.817 |        0.842 | 3/3               |
| scalp               |      60 |      0.856 |        0.85  | 3/3               |
| synth_krig          |      10 |      0.559 |        0.539 | 1/3               |
| synth_krig          |      30 |      0.589 |        0.583 | 2/3               |
| synth_krig          |      60 |      0.633 |        0.583 | 1/3               |
| synth_krig_selftest |      10 |      0.667 |        0.689 | 2/3               |
| synth_krig_selftest |      30 |      0.75  |        0.8   | 2/3               |
| synth_krig_selftest |      60 |      0.778 |        0.883 | 2/3               |
| synth_spline        |      10 |      0.5   |        0.475 | 1/3               |
| synth_spline        |      30 |      0.5   |        0.483 | 0/3               |
| synth_spline        |      60 |      0.517 |        0.533 | 1/3               |


## Paired Wilcoxon on per-subject accuracy

| comparison | window | mean diff | p |
|---|---|---|---|
| synth_krig - synth_spline | 10s | +0.059 | 0.25 |
| synth_krig - real_ear | 10s | -0.032 | 0.5 |
| synth_spline - real_ear | 10s | -0.092 | 0.25 |
| synth_krig_selftest - real_ear | 10s | +0.075 | 0.5 |
| synth_krig - chance | 10s | +0.059 | 0.5 |
| synth_krig - synth_spline | 30s | +0.089 | 0.25 |
| synth_krig - real_ear | 30s | -0.039 | 0.5 |
| synth_spline - real_ear | 30s | -0.128 | 0.25 |
| synth_krig_selftest - real_ear | 30s | +0.122 | 0.25 |
| synth_krig - chance | 30s | +0.089 | 0.5 |
| synth_krig - synth_spline | 60s | +0.117 | 0.25 |
| synth_krig - real_ear | 60s | -0.011 | 1 |
| synth_spline - real_ear | 60s | -0.128 | 0.25 |
| synth_krig_selftest - real_ear | 60s | +0.133 | 0.25 |
| synth_krig - chance | 60s | +0.133 | 0.25 |
