"""Shared numerical helpers for skinflame."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def to_array(x) -> np.ndarray:
    """Coerce list/Series/array to a 1-D numpy float array."""
    return np.asarray(x, dtype=float).ravel()


def design_matrix(df: pd.DataFrame, cols: list[str], add_intercept: bool = True) -> np.ndarray:
    """Build a 2-D float design matrix from named columns of a DataFrame."""
    if not cols:
        X = np.empty((len(df), 0))
    else:
        X = df[cols].to_numpy(dtype=float, copy=True)
    if add_intercept:
        X = np.column_stack([np.ones(len(df)), X])
    return X


def ols(y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Plain OLS via lstsq.

    Returns (beta, cov_beta, sigma2). Cov is the classical (X'X)^{-1} sigma^2.
    """
    n, p = X.shape
    XtX = X.T @ X
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    dof = max(n - p, 1)
    sigma2 = float(resid @ resid) / dof
    cov = XtX_inv * sigma2
    return beta, cov, sigma2


def percentile_ci(samples: np.ndarray, alpha: float = 0.05) -> tuple[float, float]:
    lo = float(np.quantile(samples, alpha / 2))
    hi = float(np.quantile(samples, 1 - alpha / 2))
    return lo, hi


def bca_ci(samples: np.ndarray, theta_hat: float, alpha: float = 0.05) -> tuple[float, float]:
    """Bias-corrected and accelerated bootstrap CI.

    Acceleration is approximated with the empirical influence via jackknife on
    the bootstrap samples themselves (a common pragmatic shortcut when the
    original statistic is expensive).
    """
    samples = np.asarray(samples, dtype=float)
    n = samples.size
    # Bias correction
    p0 = np.mean(samples < theta_hat)
    if p0 in (0.0, 1.0):
        return percentile_ci(samples, alpha)
    z0 = stats.norm.ppf(p0)
    # Jackknife acceleration on the bootstrap distribution (approximation)
    jk = np.array([np.mean(np.delete(samples, i)) for i in range(n)])
    jk_mean = np.mean(jk)
    num = np.sum((jk_mean - jk) ** 3)
    den = 6.0 * (np.sum((jk_mean - jk) ** 2) ** 1.5 + 1e-12)
    a = num / den
    z_lo = stats.norm.ppf(alpha / 2)
    z_hi = stats.norm.ppf(1 - alpha / 2)
    a1 = stats.norm.cdf(z0 + (z0 + z_lo) / (1 - a * (z0 + z_lo)))
    a2 = stats.norm.cdf(z0 + (z0 + z_hi) / (1 - a * (z0 + z_hi)))
    lo = float(np.quantile(samples, np.clip(a1, 1e-4, 1 - 1e-4)))
    hi = float(np.quantile(samples, np.clip(a2, 1e-4, 1 - 1e-4)))
    return lo, hi


def fdr_bh(pvals: np.ndarray, q: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg FDR.

    Returns (rejected, q_values) — q_values are BH-adjusted p-values.
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity from the bottom
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    q_values = np.empty_like(adj)
    q_values[order] = adj
    rejected = q_values <= q
    return rejected, q_values


def standardize(X: np.ndarray) -> np.ndarray:
    """Column-wise z-score; constant columns become zero."""
    X = np.asarray(X, dtype=float)
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (X - mu) / sd
