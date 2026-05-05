"""Plotting helpers (optional; only used by the demo + notebooks)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def effects_dataframe(result) -> pd.DataFrame:
    """Tidy DataFrame of effects (one row per effect) for plotting."""
    rows = []
    for label, est in [
        ("total", result.total_effect),
        ("direct", result.direct_effect),
        ("indirect", result.indirect_effect),
    ]:
        rows.append(
            {
                "effect": label,
                "estimate": est.estimate,
                "ci_low": est.ci_low,
                "ci_high": est.ci_high,
            }
        )
    for blk, est in result.per_block.items():
        rows.append(
            {
                "effect": f"indirect[{blk}]",
                "estimate": est.estimate,
                "ci_low": est.ci_low,
                "ci_high": est.ci_high,
            }
        )
    return pd.DataFrame(rows)


def forest_plot(result, ax=None, title: str | None = None):
    """Matplotlib forest plot of effects with bootstrap CIs."""
    import matplotlib.pyplot as plt

    df = effects_dataframe(result)
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 0.55 * len(df) + 1.2))
    y = np.arange(len(df))[::-1]
    ax.errorbar(
        df["estimate"],
        y,
        xerr=[df["estimate"] - df["ci_low"], df["ci_high"] - df["estimate"]],
        fmt="o",
        capsize=4,
        color="#1b4965",
    )
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["effect"])
    ax.set_xlabel("Estimate (95% CI)")
    if title:
        ax.set_title(title)
    return ax
