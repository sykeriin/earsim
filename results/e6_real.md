# E6 (real data): placement search for attention decoding, window 30s

Pairs chosen on the other subjects (leave-one-subject-out), scored on the held-out subject.

Accuracy is the fraction of correctly classified windows; chance = 0.5.

|                  |   greedy_acc |   random_acc |
|:-----------------|-------------:|-------------:|
| ('cEEGrid_L', 1) |        0.546 |        0.539 |
| ('cEEGrid_L', 2) |        0.573 |        0.561 |
| ('cEEGrid_L', 3) |        0.587 |        0.567 |
| ('cEEGrid_R', 1) |        0.57  |        0.556 |
| ('cEEGrid_R', 2) |        0.594 |        0.573 |
| ('cEEGrid_R', 3) |        0.623 |        0.582 |
| ('in_ear_L', 1)  |        0.551 |        0.519 |
| ('in_ear_L', 2)  |        0.527 |        0.518 |
| ('in_ear_L', 3)  |        0.517 |        0.529 |
| ('in_ear_R', 1)  |        0.531 |        0.511 |
| ('in_ear_R', 2)  |        0.514 |        0.516 |
| ('in_ear_R', 3)  |        0.509 |        0.519 |


| device | k | greedy - random | wins | p |
|---|---|---|---|---|
| cEEGrid_L | 1 | +0.007 | 10/15 | 0.45 |
| cEEGrid_L | 2 | +0.013 | 11/15 | 0.14 |
| cEEGrid_L | 3 | +0.020 | 11/15 | 0.022 |
| cEEGrid_R | 1 | +0.014 | 9/15 | 0.33 |
| cEEGrid_R | 2 | +0.021 | 11/15 | 0.21 |
| cEEGrid_R | 3 | +0.041 | 13/15 | 0.0084 |
| in_ear_L | 1 | +0.032 | 13/15 | 0.018 |
| in_ear_L | 2 | +0.008 | 10/15 | 0.33 |
| in_ear_L | 3 | -0.012 | 6/15 | 0.28 |
| in_ear_R | 1 | +0.020 | 12/15 | 0.03 |
| in_ear_R | 2 | -0.003 | 7/15 | 0.72 |
| in_ear_R | 3 | -0.009 | 5/15 | 0.17 |

Ranking stability (mean pairwise Spearman of leave-one-subject-out single-pair rankings):

| device    |   rank_stability_rho |
|:----------|---------------------:|
| cEEGrid_L |                 0.98 |
| cEEGrid_R |                 0.96 |
| in_ear_L  |                 0.93 |
| in_ear_R  |                 0.92 |


Top single pairs by mean leave-one-trial-out correlation (all subjects):

| device    | pair       |   mean_loto_r |
|:----------|:-----------|--------------:|
| cEEGrid_L | cEL5-cEL8  |         0.017 |
| cEEGrid_L | cEL5-cEL9  |         0.017 |
| cEEGrid_L | cEL1-cEL5  |         0.015 |
| cEEGrid_L | cEL5-cEL7  |         0.015 |
| cEEGrid_L | cEL4-cEL9  |         0.015 |
| cEEGrid_R | cER3-cER6  |         0.019 |
| cEEGrid_R | cER3-cER7  |         0.018 |
| cEEGrid_R | cER2-cER6  |         0.018 |
| cEEGrid_R | cER2-cER7  |         0.017 |
| cEEGrid_R | cER6-cER10 |         0.017 |
| in_ear_L  | ELE-ELT    |         0.007 |
| in_ear_L  | ELI-ELT    |         0.005 |
| in_ear_L  | ELC-ELI    |         0.005 |
| in_ear_L  | ELE-ELI    |         0.005 |
| in_ear_L  | ELC-ELE    |         0.004 |
| in_ear_R  | ERE-ERT    |         0.006 |
| in_ear_R  | ERA-ERT    |         0.003 |
| in_ear_R  | ERB-ERT    |         0.003 |
| in_ear_R  | ERB-ERE    |         0.003 |
| in_ear_R  | ERA-ERI    |         0.003 |
