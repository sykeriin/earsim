"""Empirical-Bayes sensor covariance model and its kriging predictor.

Model (per frequency band):  Sigma = sum_k theta_k S_k + theta_0 I,  theta >= 0,
with S_k the source-group basis covariances from `sources.py`. The theta are fitted
by maximum likelihood to the *sensor* covariances of many subjects (scalp data only),
with a free scalar gain per subject profiled out so subjects of different amplitude
pool. Prediction of any target set from any source set is the Gaussian conditional
mean  W = Sigma_ts pinv(Sigma_ss),  and sampling uses the conditional covariance.

MNE with an identity source prior is the special case of equal theta over cortical
vertices; this model strictly contains it.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def mean_zero_basis(n: int) -> np.ndarray:
    """Orthonormal basis (n x n-1) of the subspace orthogonal to the all-ones vector."""
    Q, _ = np.linalg.qr(np.c_[np.ones(n), np.eye(n)[:, : n - 1]])
    return Q[:, 1:]


def ref_operator(n: int, kept: np.ndarray) -> np.ndarray:
    """Re-reference n channels to the mean of those in boolean mask `kept`."""
    R = np.eye(n)
    R -= np.outer(np.ones(n), kept / kept.sum())
    return R


def fit_theta(S: np.ndarray, covs: np.ndarray, ridge: float = 1e-3,
              maxiter: int = 400, verbose: bool = False) -> np.ndarray:
    """ML fit of theta (K+1,) from basis S (K, m, m) and data covariances covs (N, m, m).

    All matrices must already live in the same m-dim mean-zero subspace.
    The last parameter is the white-noise variance (identity basis).
    A small ridge on log-theta (towards the mean of log-theta) keeps the
    ~100 weakly identified group variances from running to 0/inf.
    """
    K, m, _ = S.shape
    N = covs.shape[0]
    Sfull = np.concatenate([S, np.eye(m)[None]])  # (K+1, m, m)
    covs = covs / np.trace(covs, axis1=1, axis2=2)[:, None, None] * m  # scale-free
    u0 = np.log(np.full(K + 1, 1.0 / (K + 1)))

    def fun(u):
        th = np.exp(u)
        Sig = np.tensordot(th, Sfull, axes=1)
        try:
            L = np.linalg.cholesky(Sig)
        except np.linalg.LinAlgError:
            return 1e10, np.zeros_like(u)
        logdet = 2 * np.log(np.diag(L)).sum()
        A = np.linalg.inv(Sig)
        AC = np.einsum("ij,njk->nik", A, covs)           # A C_s
        q = np.einsum("nii->n", AC)                       # tr(A C_s)
        M = np.einsum("nij,jk->nik", AC, A)               # A C_s A
        f = N * logdet + m * np.log(q / m).sum()
        trAS = np.einsum("ij,kji->k", A, Sfull)           # tr(A S_k)
        SM = np.einsum("kij,nij->nk", Sfull, M)           # <S_k, A C_s A>
        g = th * (N * trAS - (m / q)[:, None].T @ SM).ravel()
        mu = u[:K].mean()
        f += 0.5 * ridge * N * np.sum((u[:K] - mu) ** 2)
        g[:K] += ridge * N * (u[:K] - mu)
        return f, g

    res = minimize(fun, u0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
    if verbose:
        print(f"  fit: f={res.fun:.2f} iters={res.nit} {res.message}")
    return np.exp(res.x)


def sigma_full(S_full: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Full-channel model covariance (reference-free) from bases over all channels."""
    n = S_full.shape[1]
    return np.tensordot(theta[:-1], S_full, axes=1) + theta[-1] * np.eye(n)


def predict_operator(Sigma: np.ndarray, src_idx: np.ndarray, tgt_idx: np.ndarray,
                     rcond: float = 1e-8) -> np.ndarray:
    """W (n_tgt x n_src): conditional-mean predictor, references = mean of source channels."""
    n = Sigma.shape[0]
    kept = np.zeros(n)
    kept[src_idx] = 1
    R = ref_operator(n, kept)
    Sr = R @ Sigma @ R.T
    Sss = Sr[np.ix_(src_idx, src_idx)]
    Sts = Sr[np.ix_(tgt_idx, src_idx)]
    return Sts @ np.linalg.pinv(Sss, rcond=rcond, hermitian=True)


def conditional_cov(Sigma: np.ndarray, src_idx, tgt_idx) -> np.ndarray:
    n = Sigma.shape[0]
    kept = np.zeros(n)
    kept[src_idx] = 1
    R = ref_operator(n, kept)
    Sr = R @ Sigma @ R.T
    Sss = Sr[np.ix_(src_idx, src_idx)]
    Sts = Sr[np.ix_(tgt_idx, src_idx)]
    Stt = Sr[np.ix_(tgt_idx, tgt_idx)]
    return Stt - Sts @ np.linalg.pinv(Sss, rcond=1e-8, hermitian=True) @ Sts.T


def metrics_from_cov(W: np.ndarray, C: np.ndarray, src_idx, tgt_idx) -> dict[str, np.ndarray]:
    """Per-target r, relative RMSE and gain of pred = W x_src, from the 64-ch data covariance.

    C is the data covariance over all channels (any common reference); it is
    re-referenced to the mean of the source channels, matching how W was built.
    """
    n = C.shape[0]
    kept = np.zeros(n)
    kept[src_idx] = 1
    R = ref_operator(n, kept)
    Cr = R @ C @ R.T
    Css = Cr[np.ix_(src_idx, src_idx)]
    Cst = Cr[np.ix_(src_idx, tgt_idx)]
    vt = np.diag(Cr)[tgt_idx]
    cov = np.einsum("ti,it->t", W, Cst)
    vp = np.einsum("ti,ij,tj->t", W, Css, W)
    r = cov / np.sqrt(vp * vt)
    return {"r": r, "rel_rmse": np.sqrt(np.maximum(vp + vt - 2 * cov, 0) / vt),
            "gain": np.sqrt(vp / vt)}
