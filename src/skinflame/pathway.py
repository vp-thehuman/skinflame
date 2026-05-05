"""Pathway-level mediation.

Aggregate per-gene (or per-metabolite) mediators into a pathway score, then
run the classic single-mediator mediation per pathway plus a multi-pathway
joint analysis. Two aggregation strategies are supported:

- ``"pca"`` (default): the first principal component of the standardized
  pathway matrix, with sign chosen so the loading is positively correlated
  with the mean.
- ``"mean"``: simple mean of column z-scores.

Pathway membership is supplied as ``pathways = {pathway_name: [gene/metab cols]}``.
Genes can appear in multiple pathways (overlap is allowed).
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from .core import MediationResult, bootstrap_mediation
from .utils import fdr_bh, standardize


def _score_pathway(values: np.ndarray, method: str = "pca") -> np.ndarray:
    """Return a 1-D pathway score from a 2-D (n x p) standardized block."""
    if values.shape[1] == 0:
        return np.zeros(values.shape[0])
    if method == "mean":
        return values.mean(axis=1)
    if method == "pca":
        pca = PCA(n_components=1)
        pc1 = pca.fit_transform(values)[:, 0]
        # Sign-align with the mean so direction is interpretable
        m = values.mean(axis=1)
        if np.corrcoef(pc1, m)[0, 1] < 0:
            pc1 = -pc1
        return pc1
    raise ValueError(f"Unknown aggregation method: {method!r}")


def _build_pathway_frame(
    data: pd.DataFrame,
    pathways: Mapping[str, Sequence[str]],
    method: str = "pca",
) -> tuple[pd.DataFrame, list[str]]:
    """Return (DataFrame with one column per pathway, list of pathway col names)."""
    out = pd.DataFrame(index=data.index.copy())
    pwy_cols: list[str] = []
    for pwy, members in pathways.items():
        members = [m for m in members if m in data.columns]
        if not members:
            continue
        block = standardize(data[members].to_numpy(dtype=float))
        score = _score_pathway(block, method=method)
        col = f"pwy::{pwy}"
        out[col] = score
        pwy_cols.append(col)
    return out, pwy_cols


def pathway_mediation(
    data: pd.DataFrame,
    exposure: str,
    outcome: str,
    pathways: Mapping[str, Sequence[str]],
    confounders: Sequence[str] = (),
    n_boot: int = 1000,
    ci: str = "percentile",
    alpha: float = 0.05,
    aggregation: str = "pca",
    fdr_q: float = 0.05,
    random_state: int | None = None,
    n_jobs: int = 1,
    **_: object,
) -> MediationResult:
    """Pathway-level mediation analysis.

    Each pathway becomes a single (aggregate) mediator. Returns a
    `MediationResult` whose ``per_block`` is keyed by pathway name and whose
    ``extras["pathway_table"]`` holds per-pathway alpha/beta/indirect estimates
    with FDR-adjusted joint p-values.
    """
    confs = list(confounders)
    pwy_df, pwy_cols = _build_pathway_frame(data, pathways, method=aggregation)
    if not pwy_cols:
        raise ValueError("No pathways had any members present in `data`.")

    cols_to_carry = [exposure, outcome, *confs]
    cols_to_carry = [c for c in cols_to_carry if c in data.columns]
    full = pd.concat([data[cols_to_carry].reset_index(drop=True), pwy_df.reset_index(drop=True)], axis=1)

    # Map to a blocks dict keyed by pathway name (without the "pwy::" prefix)
    blocks = {col.removeprefix("pwy::"): [col] for col in pwy_cols}

    result = bootstrap_mediation(
        data=full,
        exposure=exposure,
        outcome=outcome,
        mediators=blocks,
        confounders=confs,
        n_boot=n_boot,
        ci=ci,
        alpha=alpha,
        random_state=random_state,
        n_jobs=n_jobs,
    )
    result.method = "pathway"

    # Build per-pathway joint-significance table from per_mediator rows
    if result.per_mediator is not None:
        tbl = result.per_mediator.copy()
        tbl["pathway"] = tbl["mediator"].str.removeprefix("pwy::")
        rejected, qvals = fdr_bh(tbl["joint_sig_p"].to_numpy(), q=fdr_q)
        tbl["q_value"] = qvals
        tbl["fdr_significant"] = rejected
        tbl = tbl[
            [
                "pathway",
                "alpha",
                "beta_M",
                "indirect",
                "p_alpha",
                "p_beta",
                "joint_sig_p",
                "q_value",
                "fdr_significant",
            ]
        ].sort_values("joint_sig_p").reset_index(drop=True)
        result.extras["pathway_table"] = tbl
    result.extras["aggregation"] = aggregation
    result.extras["fdr_q"] = fdr_q
    return result
