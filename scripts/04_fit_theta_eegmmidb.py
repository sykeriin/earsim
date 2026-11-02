"""Fit parcel variances on EEGMMIDB (all subjects, 64 ch) for reuse on other montages.

Output data/cache/theta-eegmmidb.npz:
  v[band, K]  group variance in absolute lead-field units (theta_k / trace_k)
  noise[band] white-noise variance in the same units
  names[K], bands
Data covariances are scaled to unit per-channel variance (trace = m) per subject, so
`noise` is a *relative* level; target datasets refit one noise scale on their own scalp data.
"""
import time
import numpy as np
from earsim import cache, kriging as kr, layouts, sources

z = np.load(cache.CACHE_DIR / "cov-eegbci.npz")
C, bands, names = z["C"], [str(b) for b in z["bands"]], [str(n) for n in z["ch_names"]]
pos = layouts.standard(names)
bn, S_raw = sources.cortical_raw(pos)
t = sources.avgref_trace(S_raw)
S_unit = S_raw / t[:, None, None]
B = kr.mean_zero_basis(len(names))
S_B = np.stack([B.T @ s @ B for s in S_unit])
V, N = [], []
t0 = time.time()
for bi, band in enumerate(bands):
    Cb = np.stack([B.T @ C[s, bi] @ B for s in range(C.shape[0])])
    th = kr.fit_theta(S_B, Cb)
    V.append(th[:-1] / t)
    N.append(th[-1])
    print(f"{band}: fitted ({time.time() - t0:.0f}s) noise={th[-1]:.3g}", flush=True)
np.savez(cache.CACHE_DIR / "theta-eegmmidb.npz", v=np.array(V), noise=np.array(N),
         names=np.array(bn), bands=np.array(bands))
print("saved")
