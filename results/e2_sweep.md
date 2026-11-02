# E2: scalp -> real ear-EEG (Ear-SAAD), leave-one-subject-out, n=11 subjects

Ear electrode positions are template approximations (no coordinates ship with the dataset).
Median over subjects of mean r per output set (reference-free outputs).


## geometry r=0.025|c=0.006


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.544 |    0.521 |
| beta   |         0.304 |    0.346 |
| broad  |         0.358 |    0.35  |
| delta  |         0.363 |    0.243 |
| gamma  |         0.168 |    0.188 |
| theta  |         0.466 |    0.478 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.877 |    0.887 |
| beta   |         0.708 |    0.726 |
| broad  |         0.809 |    0.821 |
| delta  |         0.788 |    0.825 |
| gamma  |         0.513 |    0.543 |
| theta  |         0.868 |    0.88  |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.652 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.104 |    0.129 |
| beta   |         0.053 |    0.015 |
| broad  |         0.116 |    0.119 |
| delta  |         0.099 |    0.133 |
| gamma  |         0.021 |   -0.015 |
| theta  |         0.111 |    0.137 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.018 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


### Paired Wilcoxon, krig_cortex vs spline (per output set and band)

| outputs | band | median dr | wins | p |
|---|---|---|---|---|
| around_local | delta | +0.131 | 10/11 | 0.0029 |
| around_local | theta | -0.004 | 4/11 | 0.41 |
| around_local | alpha | +0.000 | 6/11 | 0.7 |
| around_local | beta | -0.020 | 1/11 | 0.019 |
| around_local | gamma | -0.035 | 1/11 | 0.019 |
| around_local | broad | -0.002 | 4/11 | 0.64 |
| around_cross | delta | -0.032 | 1/11 | 0.002 |
| around_cross | theta | -0.012 | 1/11 | 0.002 |
| around_cross | alpha | -0.005 | 4/11 | 0.15 |
| around_cross | beta | -0.013 | 0/11 | 0.00098 |
| around_cross | gamma | -0.019 | 0/11 | 0.00098 |
| around_cross | broad | -0.020 | 0/11 | 0.00098 |
| in_cross | delta | -0.040 | 0/11 | 0.00098 |
| in_cross | theta | -0.020 | 0/11 | 0.00098 |
| in_cross | alpha | -0.022 | 2/11 | 0.014 |
| in_cross | beta | -0.021 | 1/11 | 0.002 |
| in_cross | gamma | -0.016 | 1/11 | 0.0049 |
| in_cross | broad | -0.027 | 0/11 | 0.00098 |
| around_in | delta | -0.028 | 1/11 | 0.0029 |
| around_in | theta | -0.021 | 0/11 | 0.00098 |
| around_in | alpha | -0.014 | 3/11 | 0.032 |
| around_in | beta | -0.007 | 5/11 | 0.83 |
| around_in | gamma | -0.002 | 5/11 | 0.9 |
| around_in | broad | -0.015 | 2/11 | 0.054 |
| in_local | delta | -0.012 | 1/11 | 0.042 |
| in_local | theta | -0.000 | 5/11 | 0.97 |
| in_local | alpha | +0.012 | 10/11 | 0.0029 |
| in_local | beta | +0.009 | 8/11 | 0.032 |
| in_local | gamma | +0.001 | 6/11 | 0.76 |
| in_local | broad | -0.001 | 5/11 | 0.24 |

## geometry r=0.025|c=0.012


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.534 |    0.534 |
| beta   |         0.291 |    0.339 |
| broad  |         0.371 |    0.377 |
| delta  |         0.401 |    0.378 |
| gamma  |         0.159 |    0.192 |
| theta  |         0.464 |    0.485 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.879 |    0.869 |
| beta   |         0.717 |    0.729 |
| broad  |         0.797 |    0.806 |
| delta  |         0.773 |    0.772 |
| gamma  |         0.528 |    0.541 |
| theta  |         0.858 |    0.863 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.15  |    0.199 |
| beta   |         0.056 |    0.07  |
| broad  |         0.112 |    0.141 |
| delta  |         0.127 |    0.11  |
| gamma  |         0.046 |    0.018 |
| theta  |         0.168 |    0.178 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.025|c=0.02


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.495 |    0.5   |
| beta   |         0.283 |    0.33  |
| broad  |         0.378 |    0.395 |
| delta  |         0.416 |    0.422 |
| gamma  |         0.168 |    0.197 |
| theta  |         0.465 |    0.475 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.856 |    0.841 |
| beta   |         0.721 |    0.723 |
| broad  |         0.774 |    0.768 |
| delta  |         0.738 |    0.74  |
| gamma  |         0.545 |    0.514 |
| theta  |         0.83  |    0.832 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.171 |    0.219 |
| beta   |         0.093 |    0.085 |
| broad  |         0.127 |    0.118 |
| delta  |         0.096 |    0.077 |
| gamma  |         0.05  |    0.041 |
| theta  |         0.166 |    0.186 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.032|c=0.006


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.555 |    0.53  |
| beta   |         0.306 |    0.346 |
| broad  |         0.362 |    0.354 |
| delta  |         0.363 |    0.254 |
| gamma  |         0.171 |    0.187 |
| theta  |         0.47  |    0.481 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.889 |    0.885 |
| beta   |         0.708 |    0.734 |
| broad  |         0.814 |    0.825 |
| delta  |         0.798 |    0.833 |
| gamma  |         0.525 |    0.544 |
| theta  |         0.874 |    0.879 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.112 |    0.123 |
| beta   |        -0.008 |    0.026 |
| broad  |         0.082 |    0.095 |
| delta  |         0.101 |    0.126 |
| gamma  |        -0.004 |   -0.008 |
| theta  |         0.138 |    0.154 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.032|c=0.012


**around_local**

| band   |   MNE |   krig_cortex |   krig_xfer |   oracle_lsq |   spline |
|:-------|------:|--------------:|------------:|-------------:|---------:|
| alpha  | 0.493 |         0.553 |       0.509 |        0.593 |    0.521 |
| beta   | 0.307 |         0.297 |       0.288 |        0.314 |    0.341 |
| broad  | 0.304 |         0.372 |       0.358 |        0.463 |    0.375 |
| delta  | 0.272 |         0.368 |       0.393 |        0.52  |    0.349 |
| gamma  | 0.184 |         0.173 |       0.132 |        0.165 |    0.197 |
| theta  | 0.448 |         0.468 |       0.466 |        0.555 |    0.478 |


**around_cross**

| band   |   MNE |   krig_cortex |   krig_xfer |   oracle_lsq |   spline |
|:-------|------:|--------------:|------------:|-------------:|---------:|
| alpha  | 0.869 |         0.876 |       0.875 |        0.892 |    0.864 |
| beta   | 0.69  |         0.715 |       0.685 |        0.756 |    0.734 |
| broad  | 0.804 |         0.799 |       0.796 |        0.817 |    0.802 |
| delta  | 0.782 |         0.776 |       0.783 |        0.843 |    0.775 |
| gamma  | 0.541 |         0.536 |       0.465 |        0.495 |    0.538 |
| theta  | 0.869 |         0.855 |       0.866 |        0.899 |    0.858 |


**in_cross**

| band   |   MNE |   krig_cortex |   krig_xfer |   oracle_lsq |   spline |
|:-------|------:|--------------:|------------:|-------------:|---------:|
| alpha  | 0.809 |         0.777 |       0.801 |        0.795 |    0.81  |
| beta   | 0.636 |         0.646 |       0.623 |        0.69  |    0.645 |
| broad  | 0.645 |         0.619 |       0.632 |        0.672 |    0.646 |
| delta  | 0.622 |         0.576 |       0.595 |        0.664 |    0.615 |
| gamma  | 0.431 |         0.463 |       0.473 |        0.401 |    0.486 |
| theta  | 0.764 |         0.729 |       0.751 |        0.767 |    0.758 |


**around_in**

| band   |   MNE |   krig_cortex |   krig_xfer |   oracle_lsq |   spline |
|:-------|------:|--------------:|------------:|-------------:|---------:|
| alpha  | 0.176 |         0.143 |       0.153 |        0.231 |    0.184 |
| beta   | 0.051 |         0.047 |       0.052 |        0.115 |    0.064 |
| broad  | 0.126 |         0.109 |       0.106 |        0.159 |    0.116 |
| delta  | 0.13  |         0.125 |       0.124 |        0.108 |    0.122 |
| gamma  | 0     |         0.007 |       0.005 |        0.083 |    0.028 |
| theta  | 0.165 |         0.151 |       0.152 |        0.143 |    0.177 |


**in_local**

| band   |    MNE |   krig_cortex |   krig_xfer |   oracle_lsq |   spline |
|:-------|-------:|--------------:|------------:|-------------:|---------:|
| alpha  |  0.011 |         0.019 |       0.016 |        0.101 |    0.004 |
| beta   |  0.008 |         0.006 |       0.012 |        0.071 |    0.001 |
| broad  |  0.019 |         0.014 |       0.01  |        0.086 |    0.01  |
| delta  |  0.023 |         0.009 |       0.01  |        0.055 |    0.019 |
| gamma  | -0.004 |        -0     |       0.001 |        0.04  |   -0.005 |
| theta  |  0.018 |         0.011 |       0.011 |        0.065 |    0.014 |


## geometry r=0.032|c=0.02


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.522 |    0.502 |
| beta   |         0.279 |    0.33  |
| broad  |         0.378 |    0.391 |
| delta  |         0.395 |    0.397 |
| gamma  |         0.165 |    0.191 |
| theta  |         0.478 |    0.469 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.848 |    0.836 |
| beta   |         0.715 |    0.725 |
| broad  |         0.769 |    0.757 |
| delta  |         0.739 |    0.742 |
| gamma  |         0.538 |    0.498 |
| theta  |         0.826 |    0.829 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.152 |    0.207 |
| beta   |         0.087 |    0.078 |
| broad  |         0.129 |    0.115 |
| delta  |         0.103 |    0.088 |
| gamma  |         0.037 |    0.04  |
| theta  |         0.163 |    0.184 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.04|c=0.006


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.548 |    0.52  |
| beta   |         0.32  |    0.345 |
| broad  |         0.37  |    0.365 |
| delta  |         0.373 |    0.341 |
| gamma  |         0.181 |    0.19  |
| theta  |         0.475 |    0.479 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.878 |    0.872 |
| beta   |         0.708 |    0.739 |
| broad  |         0.814 |    0.818 |
| delta  |         0.794 |    0.791 |
| gamma  |         0.534 |    0.543 |
| theta  |         0.87  |    0.873 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.123 |    0.13  |
| beta   |         0.017 |    0.046 |
| broad  |         0.095 |    0.105 |
| delta  |         0.11  |    0.125 |
| gamma  |         0     |    0.017 |
| theta  |         0.143 |    0.146 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.04|c=0.012


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.556 |    0.526 |
| beta   |         0.314 |    0.337 |
| broad  |         0.37  |    0.377 |
| delta  |         0.365 |    0.355 |
| gamma  |         0.185 |    0.194 |
| theta  |         0.481 |    0.479 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.858 |    0.853 |
| beta   |         0.712 |    0.738 |
| broad  |         0.798 |    0.796 |
| delta  |         0.771 |    0.77  |
| gamma  |         0.531 |    0.527 |
| theta  |         0.852 |    0.854 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.14  |    0.155 |
| beta   |         0.036 |    0.07  |
| broad  |         0.111 |    0.12  |
| delta  |         0.131 |    0.113 |
| gamma  |         0.015 |    0.033 |
| theta  |         0.137 |    0.163 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |


## geometry r=0.04|c=0.02


**around_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.531 |    0.498 |
| beta   |         0.288 |    0.328 |
| broad  |         0.377 |    0.385 |
| delta  |         0.398 |    0.403 |
| gamma  |         0.177 |    0.191 |
| theta  |         0.479 |    0.469 |


**around_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.833 |    0.827 |
| beta   |         0.707 |    0.727 |
| broad  |         0.769 |    0.755 |
| delta  |         0.748 |    0.746 |
| gamma  |         0.513 |    0.486 |
| theta  |         0.825 |    0.824 |


**in_cross**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.777 |    0.81  |
| beta   |         0.646 |    0.645 |
| broad  |         0.619 |    0.646 |
| delta  |         0.576 |    0.615 |
| gamma  |         0.463 |    0.486 |
| theta  |         0.729 |    0.758 |


**around_in**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.145 |    0.184 |
| beta   |         0.069 |    0.076 |
| broad  |         0.121 |    0.115 |
| delta  |         0.106 |    0.1   |
| gamma  |         0.028 |    0.045 |
| theta  |         0.15  |    0.173 |


**in_local**

| band   |   krig_cortex |   spline |
|:-------|--------------:|---------:|
| alpha  |         0.019 |    0.004 |
| beta   |         0.006 |    0.001 |
| broad  |         0.014 |    0.01  |
| delta  |         0.009 |    0.019 |
| gamma  |        -0     |   -0.005 |
| theta  |         0.011 |    0.014 |

