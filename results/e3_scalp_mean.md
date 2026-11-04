# E3: zero-shot attention decoding from synthetic ear-EEG (Ear-SAAD), reference = scalp_mean

**scalp_mean reference is NOT device-realistic**: the decoder input then contains the mean scalp signal, so this row measures scalp information, not ear information.

Accuracy = fraction of windows where the attended stream is picked; chance = 0.5.

Per-subject accuracy, mean over subjects (and number of subjects above chance by binomial test p<0.05).

| method              |   win_s |   mean_acc |   median_acc | subjects>chance   |
|:--------------------|--------:|-----------:|-------------:|:------------------|
| real_ear            |      10 |      0.579 |        0.581 | 12/15             |
| real_ear            |      30 |      0.633 |        0.633 | 12/15             |
| real_ear            |      60 |      0.671 |        0.667 | 12/15             |
| scalp               |      10 |      0.676 |        0.681 | 15/15             |
| scalp               |      30 |      0.783 |        0.767 | 15/15             |
| scalp               |      60 |      0.85  |        0.85  | 15/15             |
| synth_krig          |      10 |      0.517 |        0.508 | 3/15              |
| synth_krig          |      30 |      0.535 |        0.517 | 3/15              |
| synth_krig          |      60 |      0.538 |        0.517 | 4/15              |
| synth_krig_selftest |      10 |      0.644 |        0.642 | 14/15             |
| synth_krig_selftest |      30 |      0.719 |        0.717 | 13/15             |
| synth_krig_selftest |      60 |      0.779 |        0.8   | 12/15             |
| synth_oracle_other  |      10 |      0.524 |        0.519 | 3/15              |
| synth_oracle_other  |      30 |      0.543 |        0.525 | 3/15              |
| synth_oracle_other  |      60 |      0.568 |        0.55  | 4/15              |
| synth_oracle_self   |      10 |      0.536 |        0.528 | 5/15              |
| synth_oracle_self   |      30 |      0.562 |        0.558 | 5/15              |
| synth_oracle_self   |      60 |      0.576 |        0.583 | 5/15              |
| synth_spline        |      10 |      0.505 |        0.514 | 2/15              |
| synth_spline        |      30 |      0.512 |        0.525 | 1/15              |
| synth_spline        |      60 |      0.53  |        0.533 | 3/15              |


## Paired Wilcoxon on per-subject accuracy

| comparison | window | mean diff | p |
|---|---|---|---|
| synth_krig - synth_spline | 10s | +0.012 | 0.211 |
| synth_krig - real_ear | 10s | -0.061 | 0.00236 |
| synth_spline - real_ear | 10s | -0.074 | 0.000803 |
| synth_krig_selftest - real_ear | 10s | +0.065 | 0.00286 |
| synth_oracle_other - synth_krig | 10s | +0.006 | 0.426 |
| synth_oracle_self - synth_krig | 10s | +0.019 | 0.173 |
| synth_oracle_self - real_ear | 10s | -0.043 | 0.00237 |
| synth_krig - chance | 10s | +0.017 | 0.201 |
| synth_krig - synth_spline | 30s | +0.023 | 0.169 |
| synth_krig - real_ear | 30s | -0.098 | 0.00195 |
| synth_spline - real_ear | 30s | -0.121 | 0.000985 |
| synth_krig_selftest - real_ear | 30s | +0.087 | 0.00258 |
| synth_oracle_other - synth_krig | 30s | +0.008 | 0.378 |
| synth_oracle_self - synth_krig | 30s | +0.027 | 0.28 |
| synth_oracle_self - real_ear | 30s | -0.071 | 0.0115 |
| synth_krig - chance | 30s | +0.035 | 0.0839 |
| synth_krig - synth_spline | 60s | +0.008 | 0.807 |
| synth_krig - real_ear | 60s | -0.133 | 0.00423 |
| synth_spline - real_ear | 60s | -0.141 | 0.000653 |
| synth_krig_selftest - real_ear | 60s | +0.108 | 0.0064 |
| synth_oracle_other - synth_krig | 60s | +0.030 | 0.107 |
| synth_oracle_self - synth_krig | 60s | +0.038 | 0.22 |
| synth_oracle_self - real_ear | 60s | -0.096 | 0.00628 |
| synth_krig - chance | 60s | +0.038 | 0.209 |
