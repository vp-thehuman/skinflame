"""HIMA: high-dimensional mediation analysis.

Implementation broadly follows Zhang et al. (2016) "Estimating and Testing
High-dimensional Mediation Effects in Epigenetic Studies" (Bioinformatics):

  1. Sure-Independence Screening (SIS): rank mediators by |corr(M_j, Y | X)|
     and keep the top d = floor(2 n / log n).
  2. Penalized regression of Y on the screened mediators (adjusted for X and
     C) with an L1 / MCP-style penalty. We use sklearn's `LassoCV` as a
     pragmatic stand-in for MCP (the cited HIMA-MCP requires `ncvreg`); the
     joint-significance step that follows is what controls FDR, not the
     specific penalty.
  3. Joint-significance test on the surviving mediators: for each j compute
     p_alpha (X -> M_j) and p_beta (M_j -> Y | X, others) and report
     max(p_alpha, p_beta). FDR is controlled via Benjamini-Hochberg.
  4. Refit the outcome and mediator models on the *selected* mediators only
     and run the standard classic-mediation bootstrap to obtain CIs on the
     total / direct / indirect effects.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LassoCV

from .core import MediationResult, bootstrap_mediation
from .utils import design_matrix, fdr_bh, ols, standardize


def _sis_screen(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: list[str],
    confounders: list[str],
    d: int | None = None,
) -> list[str]:
    """Sure-Independence Screening: marginal |corr(M_j, Y | X, C)|."""
    n = len(data)
    if d is None:
        d = max(int(np.floor(2 * n / max(np.log(max(n, 2)), 1.0))), 5)
    d = min(d, len(mediators))

    # Residualize Y and each M on (1, X, C) to approximate partial correlation.
    Z = design_matrix(data, [exposure, *confounders], add_intercept=True)
    y = data[outcome].to_numpy(dtype=float)
    PZ = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T  # projection
    yr = y - PZ @ y

    scores = np.zeros(len(mediators))
    for j, m in enumerate(mediators):
        mj = data[m].to_numpy(dtype=float)
        mr = mj - PZ @ mj
        denom = float(np.linalg.norm(mr) * np.linalg.norm(yr))
        scores[j] = abs(float(mr @ yr) / denom) if denom > 0 else 0.0

    order = np.argsort(-scores)
    keep = [mediators[i] for i in order[:d]]
    return keep


def _penalized_select(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: list[str],
    confounders: list[str],
    random_state: int | None = None,
) -> list[str]:
    """LassoCV on Y ~ X + M (screened) + C, return mediators with non-zero coef."""
    if not mediators:
        return []
    X_block = data[mediators].to_numpy(dtype=float)
    X_block = standardize(X_block)
    extras = data[[exposure, *confounders]].to_numpy(dtype=float)
    extras = standardize(extras) if extras.shape[1] > 0 else extras
    y = data[outcome].to_numpy(dtype=float)
    n_features_extras = extras.shape[1]

    # Fit Lasso but force exposure/confounders via two-stage residualization:
    # residualize y on extras first, then Lasso on M only.
    if n_features_extras > 0:
        Pe = extras @ np.linalg.pinv(extras.T @ extras) @ extras.T
        y_res = y - Pe @ y
        M_res = X_block - Pe @ X_block
    else:
        y_res = y
        M_res = X_block

    try:
        lasso = LassoCV(cv=min(5, max(2, len(data) // 20)), random_state=random_state, n_jobs=1)
        lasso.fit(M_res, y_res)
    except Exception:
        return mediators  # degrade gracefully
    nz = np.where(np.abs(lasso.coef_) > 1e-8)[0]
    return [mediators[i] for i in nz]


def _joint_significance(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: list[str],
    confounders: list[str],
) -> pd.DataFrame:
    """Compute p_alpha (X→M_j) and p_beta (M_j→Y | X, others) per mediator."""
    Xc = design_matrix(data, [exposure, *confounders], add_intercept=True)
    p_alpha = np.zeros(len(mediators))
    alphas = np.zeros(len(mediators))

    for j, m in enumerate(mediators):
        y = data[m].to_numpy(dtype=float)
        beta, cov, _ = ols(y, Xc)
        a = beta[1]
        se = float(np.sqrt(max(cov[1, 1], 0.0)))
        z = a / se if se > 0 else 0.0
        p_alpha[j] = 2 * (1 - stats.norm.cdf(abs(z)))
        alphas[j] = a

    # Outcome model with all selected mediators
    cols = [exposure, *mediators, *confounders]
    Xo = design_matrix(data, cols, add_intercept=True)
    yo = data[outcome].to_numpy(dtype=float)
    beta_o, cov_o, _ = ols(yo, Xo)
    n_med = len(mediators)
    bm = beta_o[2 : 2 + n_med]
    se_m = np.sqrt(np.maximum(np.diag(cov_o)[2 : 2 + n_med], 0.0))
    z_m = np.where(se_m > 0, bm / se_m, 0.0)
    p_beta = 2 * (1 - stats.norm.cdf(np.abs(z_m)))

    js = np.maximum(p_alpha, p_beta)
    return pd.DataFrame(
        {
            "mediator": mediators,
            "alpha": alphas,
            "beta_M": bm,
            "p_alpha": p_alpha,
            "p_beta": p_beta,
            "joint_sig_p": js,
            "indirect": alphas * bm,
        }
    )


def hima(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: Sequence[str] | Mapping[str, Sequence[str]],
    confounders: Sequence[str] = (),
    n_boot: int = 1000,
    ci: str = "percentile",
    alpha: float = 0.05,
    fdr_q: float = 0.05,
    sis_size: int | None = None,
    random_state: int | None = None,
    n_jobs: int = 1,
    **_: object,
) -> MediationResult:
    """High-dimensional mediation pipeline. Returns a `MediationResult` with
    additional `extras["selected_mediators"]` (DataFrame of selected mediators
    with FDR-adjusted joint p-values)."""
    if isinstance(mediators, Mapping):
        flat = [c for cols in mediators.values() for c in cols]
        blocks = {k: list(v) for k, v in mediators.items()}
    else:
        flat = list(mediators)
        blocks = {"all": flat}
    confs = list(confounders)

    # 1. Screen
    screened = _sis_screen(data, exposure, outcome, flat, confs, d=sis_size)
    # 2. Penalized selection
    selected = _penalized_select(data, exposure, outcome, screened, confs, random_state)
    if not selected:
        # fall back to top-5 from screening so downstream steps don't explode
        selected = screened[: min(5, len(screened))]

    # 3. Joint significance + FDR
    js_df = _joint_significance(data, exposure, outcome, selected, confs)
    rejected, qvals = fdr_bh(js_df["joint_sig_p"].to_numpy(), q=fdr_q)
    js_df["q_value"] = qvals
    js_df["fdr_significant"] = rejected
    js_df = js_df.sort_values("joint_sig_p").reset_index(drop=True)

    # 4. Bootstrap classic mediation on the *selected* set, preserving block
    # membership where possible.
    selected_set = set(selected)
    blocks_selected: dict[str, list[str]] = {}
    for blk, cols in blocks.items():
        keep = [c for c in cols if c in selected_set]
        if keep:
            blocks_selected[blk] = keep
    if not blocks_selected:
        blocks_selected = {"selected": selected}

    result = bootstrap_mediation(
        data=data,
        exposure=exposure,
        outcome=outcome,
        mediators=blocks_selected,
        confounders=confs,
        n_boot=n_boot,
        ci=ci,
        alpha=alpha,
        random_state=random_state,
        n_jobs=n_jobs,
    )
    result.method = "hima"
    result.extras["selected_mediators"] = js_df
    result.extras["screened_mediators"] = screened
    result.extras["fdr_q"] = fdr_q
    return result
