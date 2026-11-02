import numpy as np

from earsim import kriging as kr


def _bases(n=20, K=6, seed=0):
    rng = np.random.default_rng(seed)
    B = kr.mean_zero_basis(n)
    S = []
    for _ in range(K):
        a = rng.standard_normal((n, 3))
        S.append(B.T @ (a @ a.T) @ B)
    S = np.stack(S)
    return S / np.trace(S, axis1=1, axis2=2)[:, None, None], B


def test_mean_zero_basis():
    B = kr.mean_zero_basis(10)
    np.testing.assert_allclose(B.T @ B, np.eye(9), atol=1e-12)
    np.testing.assert_allclose(B.T @ np.ones(10), 0, atol=1e-12)


def test_fit_recovers_theta_shape():
    S, _ = _bases()
    K, m, _ = S.shape
    th_true = np.array([4.0, 0.2, 2.0, 0.1, 1.0, 3.0, 0.05])
    Sig = np.tensordot(th_true[:-1], S, 1) + th_true[-1] * np.eye(m)
    rng = np.random.default_rng(1)
    L = np.linalg.cholesky(Sig)
    covs = []
    for g in rng.uniform(0.5, 2.0, 30):                 # subjects with different gain
        X = L @ rng.standard_normal((m, 20000))
        covs.append(g * X @ X.T / X.shape[1])
    covs = np.stack(covs)
    th = kr.fit_theta(S, covs, ridge=0.0, maxiter=2000)
    # scale is not identified (per-subject gain profiled out): compare shapes
    a, b = th / th.sum(), th_true / th_true.sum()
    np.testing.assert_allclose(a, b, atol=0.03)


def test_conditional_cov_psd_and_prediction_shape():
    rng = np.random.default_rng(2)
    a = rng.standard_normal((12, 5))
    Sig = a @ a.T + 0.1 * np.eye(12)
    src, tgt = np.arange(8), np.arange(8, 12)
    W = kr.predict_operator(Sig, src, tgt)
    assert W.shape == (4, 8)
    ev = np.linalg.eigvalsh(kr.conditional_cov(Sig, src, tgt))
    assert ev.min() > -1e-9


def test_metrics_perfect_predictor():
    rng = np.random.default_rng(3)
    a = rng.standard_normal((10, 10))
    C = a @ a.T
    src, tgt = np.arange(7), np.arange(7, 10)
    # oracle regression from the data covariance itself
    W = kr.predict_operator(C, src, tgt)
    m = kr.metrics_from_cov(W, C, src, tgt)
    assert np.all(m["r"] > 0) and np.all(m["r"] <= 1 + 1e-9)
