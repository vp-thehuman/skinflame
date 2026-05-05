"""Core mediation analysis: total / direct / indirect effects with bootstrap CIs.

The estimators here follow the potential-outcomes framing of Imai et al. (2010)
for the linear-Gaussian case, where natural direct and natural indirect effects
collapse to the familiar Baron-Kenny product-of-coefficients decomposition:

    M_i = alpha0 + alpha * X_i + gamma_M' C_i + e_M
    Y_i = beta0  + beta_X * X_i + beta_M' M_i + gamma_Y' C_i + e_Y

    Total effect (TE)            = alpha * beta_M (sum) + beta_X
    Natural Indirect Effect (NIE) = sum_j alpha_j * beta_M,j
    Natural Direct Effect (NDE)   = beta_X

Single-mediator case is the special case with one alpha and one beta_M.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .utils import bca_ci, design_matrix, ols, percentile_ci


@dataclass
class EffectEstimate:
    """A single effect estimate with point estimate, SE, CI, and p-value."""
    estimate: float
    se: float | None = None
    ci_low: float | None = None
    ci_high: float | None = None
    p_value: float | None = None

    def as_row(self, name: str) -> dict:
        return {
            "effect": name,
            "estimate": self.estimate,
            "se": self.se,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "p_value": self.p_value,
        }


@dataclass
class MediationResult:
    """Container for the output of `mediation_analysis`."""

    total_effect: EffectEstimate
    direct_effect: EffectEstimate
    indirect_effect: EffectEstimate
    proportion_mediated: float | None
    per_block: dict[str, EffectEstimate] = field(default_factory=dict)
    per_mediator: pd.DataFrame | None = None
    n_boot: int = 0
    method: str = "classic"
    sobel_p: float | None = None
    joint_significance_p: float | None = None
    extras: dict = field(default_factory=dict)

    def summary(self) -> pd.DataFrame:
        rows = [
            self.total_effect.as_row("total"),
            self.direct_effect.as_row("direct"),
            self.indirect_effect.as_row("indirect"),
        ]
        for k, v in self.per_block.items():
            rows.append(v.as_row(f"indirect[{k}]"))
        df = pd.DataFrame(rows)
        df.attrs["proportion_mediated"] = self.proportion_mediated
        df.attrs["sobel_p"] = self.sobel_p
        df.attrs["joint_significance_p"] = self.joint_significance_p
        df.attrs["method"] = self.method
        df.attrs["n_boot"] = self.n_boot
        return df

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        df = self.summary()
        head = (
            f"MediationResult(method={self.method}, "
            f"prop_mediated={self.proportion_mediated:.3f}"
            if self.proportion_mediated is not None
            else f"MediationResult(method={self.method}"
        )
        return f"{head})\n{df.to_string(index=False)}"


# ---------------------------------------------------------------------------
# Point estimators
# ---------------------------------------------------------------------------

def _flatten_mediators(
    mediators: Sequence[str] | Mapping[str, Sequence[str]],
) -> tuple[list[str], dict[str, list[str]]]:
    """Return (flat_list_of_mediator_columns, blocks_dict)."""
    if isinstance(mediators, Mapping):
        blocks = {k: list(v) for k, v in mediators.items()}
        flat = [c for cols in blocks.values() for c in cols]
    else:
        flat = list(mediators)
        blocks = {"all": flat}
    # de-dup while preserving order
    seen = set()
    flat_unique = [c for c in flat if not (c in seen or seen.add(c))]
    return flat_unique, blocks


def _fit_mediator_models(
    data: pd.DataFrame,
    exposure: str,
    mediators: list[str],
    confounders: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Fit M_j = a0 + a*X + gamma'C for each mediator. Returns (alphas, alpha_se)."""
    Xc = design_matrix(data, [exposure, *confounders], add_intercept=True)
    alphas = np.zeros(len(mediators))
    alpha_se = np.zeros(len(mediators))
    for j, m in enumerate(mediators):
        y = data[m].to_numpy(dtype=float)
        beta, cov, _ = ols(y, Xc)
        alphas[j] = beta[1]  # coef on exposure
        alpha_se[j] = float(np.sqrt(max(cov[1, 1], 0.0)))
    return alphas, alpha_se


def _fit_outcome_model(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: list[str],
    confounders: list[str],
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Fit Y = b0 + b_X*X + b_M'M + gamma'C. Return (b_X, se_X, b_M, se_M)."""
    cols = [exposure, *mediators, *confounders]
    X = design_matrix(data, cols, add_intercept=True)
    y = data[outcome].to_numpy(dtype=float)
    beta, cov, _ = ols(y, X)
    bx = float(beta[1])
    se_x = float(np.sqrt(max(cov[1, 1], 0.0)))
    n_med = len(mediators)
    bm = beta[2 : 2 + n_med].astype(float)
    se_m = np.sqrt(np.maximum(np.diag(cov)[2 : 2 + n_med], 0.0))
    return bx, se_x, bm, se_m


def _fit_total_effect(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    confounders: list[str],
) -> tuple[float, float]:
    cols = [exposure, *confounders]
    X = design_matrix(data, cols, add_intercept=True)
    y = data[outcome].to_numpy(dtype=float)
    beta, cov, _ = ols(y, X)
    return float(beta[1]), float(np.sqrt(max(cov[1, 1], 0.0)))


def _point_estimates(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: list[str],
    confounders: list[str],
    blocks: dict[str, list[str]],
) -> dict:
    alphas, alpha_se = _fit_mediator_models(data, exposure, mediators, confounders)
    bx, se_x, bm, se_m = _fit_outcome_model(data, exposure, outcome, mediators, confounders)
    te, te_se = _fit_total_effect(data, exposure, outcome, confounders)

    indirect_per_med = alphas * bm
    nie = float(indirect_per_med.sum())
    nde = float(bx)
    per_block = {}
    name_to_idx = {m: i for i, m in enumerate(mediators)}
    for blk, cols in blocks.items():
        idx = [name_to_idx[c] for c in cols if c in name_to_idx]
        per_block[blk] = float(indirect_per_med[idx].sum()) if idx else 0.0

    # Sobel for the *aggregate* indirect effect: delta-method on sum_j a_j*b_j
    # Var(sum a_j b_j) approx sum (b_j^2 var(a_j) + a_j^2 var(b_j))
    var_nie = float(np.sum(bm**2 * alpha_se**2 + alphas**2 * se_m**2))
    se_nie = float(np.sqrt(max(var_nie, 0.0)))
    sobel_z = nie / se_nie if se_nie > 0 else np.nan
    sobel_p = float(2 * (1 - stats.norm.cdf(abs(sobel_z)))) if np.isfinite(sobel_z) else np.nan

    # Joint significance: for each mediator, max(p_alpha, p_beta), then take min across
    # mediators to test "any path active". For single mediator this is the standard JS test.
    z_alpha = np.where(alpha_se > 0, alphas / alpha_se, 0.0)
    z_beta = np.where(se_m > 0, bm / se_m, 0.0)
    p_alpha = 2 * (1 - stats.norm.cdf(np.abs(z_alpha)))
    p_beta = 2 * (1 - stats.norm.cdf(np.abs(z_beta)))
    js_per = np.maximum(p_alpha, p_beta)
    js_p = float(js_per.min()) if js_per.size else np.nan

    return {
        "alphas": alphas,
        "alpha_se": alpha_se,
        "bm": bm,
        "se_m": se_m,
        "bx": bx,
        "se_x": se_x,
        "te": te,
        "te_se": te_se,
        "nie": nie,
        "nde": nde,
        "per_block_nie": per_block,
        "sobel_p": sobel_p,
        "se_nie": se_nie,
        "indirect_per_med": indirect_per_med,
        "joint_significance_p": js_p,
        "p_alpha": p_alpha,
        "p_beta": p_beta,
    }


# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------

def bootstrap_mediation(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    mediators: Sequence[str] | Mapping[str, Sequence[str]],
    confounders: Sequence[str] = (),
    n_boot: int = 1000,
    ci: str = "percentile",
    alpha: float = 0.05,
    random_state: int | None = None,
    n_jobs: int = 1,
) -> MediationResult:
    """Bootstrap CIs for total / direct / indirect (per-block) effects.

    Parameters
    ----------
    data : DataFrame containing all referenced columns.
    exposure, outcome : column names.
    mediators : list of column names OR dict mapping block name -> list of cols.
    confounders : optional list of confounder column names.
    n_boot : bootstrap replicates.
    ci : "percentile" or "bca".
    alpha : two-sided alpha level (CI is 1 - alpha).
    """
    flat, blocks = _flatten_mediators(mediators)
    confs = list(confounders)
    point = _point_estimates(data, exposure, outcome, flat, confs, blocks)

    rng = np.random.default_rng(random_state)
    n = len(data)

    boots_te = np.empty(n_boot)
    boots_nde = np.empty(n_boot)
    boots_nie = np.empty(n_boot)
    boots_per_block = {k: np.empty(n_boot) for k in blocks}

    if n_jobs and n_jobs > 1:
        try:
            from joblib import Parallel, delayed
        except ImportError:  # pragma: no cover
            n_jobs = 1

    def _one(seed: int):
        rs = np.random.default_rng(seed)
        idx = rs.integers(0, n, size=n)
        sample = data.iloc[idx].reset_index(drop=True)
        try:
            est = _point_estimates(sample, exposure, outcome, flat, confs, blocks)
        except np.linalg.LinAlgError:
            return None
        return est

    if n_jobs and n_jobs > 1:
        seeds = rng.integers(0, 2**31 - 1, size=n_boot)
        results = Parallel(n_jobs=n_jobs)(delayed(_one)(int(s)) for s in seeds)
    else:
        results = [_one(int(rng.integers(0, 2**31 - 1))) for _ in range(n_boot)]

    valid = [r for r in results if r is not None]
    if len(valid) < n_boot * 0.9:
        warnings.warn(
            f"Only {len(valid)}/{n_boot} bootstrap fits succeeded; CIs may be unreliable.",
            stacklevel=2,
        )
    n_boot = len(valid)
    if n_boot == 0:
        raise RuntimeError("All bootstrap fits failed; check input data for collinearity.")

    boots_te = np.array([r["te"] for r in valid])
    boots_nde = np.array([r["nde"] for r in valid])
    boots_nie = np.array([r["nie"] for r in valid])
    boots_per_block = {
        k: np.array([r["per_block_nie"][k] for r in valid]) for k in blocks
    }

    def _wrap(theta_hat: float, samples: np.ndarray) -> EffectEstimate:
        if ci == "bca":
            lo, hi = bca_ci(samples, theta_hat, alpha)
        else:
            lo, hi = percentile_ci(samples, alpha)
        # bootstrap p-value: 2 * min(P(boot >=0), P(boot<=0)) under H0=0
        p = 2.0 * min((samples >= 0).mean(), (samples <= 0).mean())
        p = float(min(max(p, 1.0 / (len(samples) + 1)), 1.0))
        return EffectEstimate(
            estimate=float(theta_hat),
            se=float(np.std(samples, ddof=1)),
            ci_low=float(lo),
            ci_high=float(hi),
            p_value=p,
        )

    te_est = _wrap(point["te"], boots_te)
    nde_est = _wrap(point["nde"], boots_nde)
    nie_est = _wrap(point["nie"], boots_nie)
    per_block_est = {k: _wrap(point["per_block_nie"][k], v) for k, v in boots_per_block.items()}

    pm = float(point["nie"] / point["te"]) if abs(point["te"]) > 1e-9 else None

    per_med_df = pd.DataFrame(
        {
            "mediator": flat,
            "alpha": point["alphas"],
            "alpha_se": point["alpha_se"],
            "beta_M": point["bm"],
            "beta_M_se": point["se_m"],
            "indirect": point["indirect_per_med"],
            "p_alpha": point["p_alpha"],
            "p_beta": point["p_beta"],
            "joint_sig_p": np.maximum(point["p_alpha"], point["p_beta"]),
        }
    )

    return MediationResult(
        total_effect=te_est,
        direct_effect=nde_est,
        indirect_effect=nie_est,
        proportion_mediated=pm,
        per_block=per_block_est,
        per_mediator=per_med_df,
        n_boot=n_boot,
        method="classic",
        sobel_p=point["sobel_p"],
        joint_significance_p=point["joint_significance_p"],
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def mediation_analysis(
    exposure: str,
    outcome: str,
    mediators: Sequence[str] | Mapping[str, Sequence[str]],
    data: pd.DataFrame,
    confounders: Sequence[str] = (),
    method: str = "classic",
    n_boot: int = 1000,
    ci: str = "percentile",
    alpha: float = 0.05,
    dag=None,
    pathways: Mapping[str, Sequence[str]] | None = None,
    random_state: int | None = None,
    n_jobs: int = 1,
    **kwargs,
) -> MediationResult:
    """High-level mediation entry point.

    Parameters
    ----------
    method : {"classic", "hima", "pathway"}
        - "classic": OLS-based decomposition with bootstrap CIs.
        - "hima":    high-dimensional mediation; selects active mediators
                     before computing per-block / global effects. Suitable for
                     transcriptome-scale mediator blocks.
        - "pathway": aggregate mediators by pathway, then run classic on
                     pathway scores. Requires ``pathways=`` mapping.
    dag : optional networkx.DiGraph
        Graph over (exposure, outcome, mediators, confounders). If supplied,
        a back-door adjustment set is derived and used in place of/added to
        the user-supplied confounders. See :func:`skinflame.dag.adjustment_set`.
    pathways : mapping pathway -> list of mediator columns (for method="pathway").
    """
    from .dag import adjustment_set, validate_dag

    confs = list(confounders)
    if dag is not None:
        validate_dag(dag, exposure=exposure, outcome=outcome)
        adj = adjustment_set(dag, exposure=exposure, outcome=outcome)
        # Union with user confounders (preserve order)
        for c in adj:
            if c in data.columns and c != exposure and c != outcome and c not in confs:
                confs.append(c)

    if method == "classic":
        return bootstrap_mediation(
            data=data,
            exposure=exposure,
            outcome=outcome,
            mediators=mediators,
            confounders=confs,
            n_boot=n_boot,
            ci=ci,
            alpha=alpha,
            random_state=random_state,
            n_jobs=n_jobs,
        )
    if method == "hima":
        from .hima import hima

        return hima(
            data=data,
            exposure=exposure,
            outcome=outcome,
            mediators=mediators,
            confounders=confs,
            n_boot=n_boot,
            ci=ci,
            alpha=alpha,
            random_state=random_state,
            n_jobs=n_jobs,
            **kwargs,
        )
    if method == "pathway":
        from .pathway import pathway_mediation

        if pathways is None:
            raise ValueError("method='pathway' requires `pathways=` mapping.")
        return pathway_mediation(
            data=data,
            exposure=exposure,
            outcome=outcome,
            pathways=pathways,
            confounders=confs,
            n_boot=n_boot,
            ci=ci,
            alpha=alpha,
            random_state=random_state,
            n_jobs=n_jobs,
            **kwargs,
        )
    raise ValueError(f"Unknown method: {method!r}")
