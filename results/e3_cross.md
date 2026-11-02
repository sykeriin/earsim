# E3: zero-shot attention decoding from synthetic ear-EEG (Ear-SAAD), reference = cross


Accuracy = fraction of windows where the attended stream is picked; chance = 0.5.

Per-subject accuracy, mean over subjects (and number of subjects above chance by binomial test p<0.05).

| method              |   win_s |   mean_acc |   median_acc | subjects>chance   |
|:--------------------|--------:|-----------:|-------------:|:------------------|
| real_ear            |      10 |      0.579 |        0.589 | 13/15             |
| real_ear            |      30 |      0.622 |        0.617 | 11/15             |
| real_ear            |      60 |      0.679 |        0.683 | 12/15             |
| scalp               |      10 |      0.676 |        0.681 | 15/15             |
| scalp               |      30 |      0.783 |        0.767 | 15/15             |
| scalp               |      60 |      0.85  |        0.85  | 15/15             |
| synth_krig          |      10 |      0.511 |        0.506 | 2/15              |
| synth_krig          |      30 |      0.512 |        0.508 | 2/15              |
| synth_krig          |      60 |      0.51  |        0.483 | 2/15              |
| synth_krig_selftest |      10 |      0.639 |        0.644 | 12/15             |
| synth_krig_selftest |      30 |      0.713 |        0.7   | 12/15             |
| synth_krig_selftest |      60 |      0.788 |        0.8   | 12/15             |
| synth_spline        |      10 |      0.502 |        0.492 | 2/15              |
| synth_spline        |      30 |      0.497 |        0.492 | 1/15              |
| synth_spline        |      60 |      0.5   |        0.467 | 3/15              |


## Paired Wilcoxon on per-subject accuracy

| comparison | window | mean diff | p |
|---|---|---|---|
| synth_krig - synth_spline | 10s | +0.009 | 0.211 |
| synth_krig - real_ear | 10s | -0.067 | 0.00169 |
| synth_spline - real_ear | 10s | -0.076 | 0.000648 |
| synth_krig_selftest - real_ear | 10s | +0.060 | 0.00918 |
| synth_krig - chance | 10s | +0.011 | 0.362 |
| synth_krig - synth_spline | 30s | +0.016 | 0.176 |
| synth_krig - real_ear | 30s | -0.110 | 0.00122 |
| synth_spline - real_ear | 30s | -0.126 | 0.000653 |
| synth_krig_selftest - real_ear | 30s | +0.091 | 0.00757 |
| synth_krig - chance | 30s | +0.012 | 0.752 |
| synth_krig - synth_spline | 60s | +0.010 | 0.637 |
| synth_krig - real_ear | 60s | -0.169 | 0.000801 |
| synth_spline - real_ear | 60s | -0.179 | 0.00122 |
| synth_krig_selftest - real_ear | 60s | +0.109 | 0.0105 |
| synth_krig - chance | 60s | +0.010 | 0.69 |
