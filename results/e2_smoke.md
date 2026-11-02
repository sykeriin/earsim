# E2: scalp -> real ear-EEG (Ear-SAAD), leave-one-subject-out, n=3 subjects

Ear electrode positions are template approximations (no coordinates ship with the dataset).
Median over subjects of mean r per output set (reference-free outputs).


## geometry r=0.032|c=0.012


**around_local**

| band   |   MNE |   krig_cortex |   oracle_lsq |   spline |
|:-------|------:|--------------:|-------------:|---------:|
| alpha  | 0.564 |          0.57 |        0.593 |    0.547 |
| broad  | 0.255 |          0.35 |        0.481 |    0.363 |


**around_cross**

| band   |   MNE |   krig_cortex |   oracle_lsq |   spline |
|:-------|------:|--------------:|-------------:|---------:|
| alpha  | 0.869 |         0.877 |        0.916 |    0.864 |
| broad  | 0.826 |         0.821 |        0.823 |    0.825 |


**in_cross**

| band   |   MNE |   krig_cortex |   oracle_lsq |   spline |
|:-------|------:|--------------:|-------------:|---------:|
| alpha  | 0.854 |         0.841 |        0.694 |    0.85  |
| broad  | 0.645 |         0.62  |        0.642 |    0.646 |


**around_in**

| band   |   MNE |   krig_cortex |   oracle_lsq |   spline |
|:-------|------:|--------------:|-------------:|---------:|
| alpha  | 0.31  |         0.31  |        0.081 |    0.307 |
| broad  | 0.183 |         0.185 |        0.188 |    0.172 |


**in_local**

| band   |   MNE |   krig_cortex |   oracle_lsq |   spline |
|:-------|------:|--------------:|-------------:|---------:|
| alpha  | 0.011 |         0.044 |        0     |    0.023 |
| broad  | 0.073 |         0.053 |        0.154 |    0.06  |


### Paired Wilcoxon, krig_cortex vs spline (per output set and band)

| outputs | band | median dr | wins | p |
|---|---|---|---|---|
| around_local | alpha | +0.015 | 2/3 | nan |
| around_local | broad | -0.013 | 1/3 | nan |
| around_cross | alpha | +0.011 | 2/3 | nan |
| around_cross | broad | -0.004 | 1/3 | nan |
| in_cross | alpha | -0.009 | 1/3 | nan |
| in_cross | broad | -0.026 | 1/3 | nan |
| around_in | alpha | +0.003 | 2/3 | nan |
| around_in | broad | +0.012 | 3/3 | nan |
| in_local | alpha | +0.021 | 2/3 | nan |
| in_local | broad | -0.007 | 1/3 | nan |
