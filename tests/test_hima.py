import numpy as np
import pytest

from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort
from skinflame.hima import hima


def test_hima_selects_some_active_mediators():
    b = simulate_ad_cohort(n=400, p_transcript=120, p_metab=30, p_microbiome=15, seed=10)
    res = hima(
        data=b.df,
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators=b.transcript_cols + b.metab_cols + b.microbiome_cols,
        confounders=["age", "sex", "BMI"],
        n_boot=80,
        random_state=0,
    )
    sel = res.extras["selected_mediators"]
    # Top selected should include at least one of the truly active mediators
    actives = set(b.transcript_cols[:6]) | set(b.metab_cols[:3]) | set(b.microbiome_cols[:2])
    top_picks = set(sel.head(10)["mediator"].tolist())
    assert len(top_picks & actives) >= 1


def test_hima_via_high_level_api():
    b = simulate_ad_cohort(n=300, p_transcript=80, p_metab=20, p_microbiome=10, seed=11)
    res = mediation_analysis(
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators={"tx": b.transcript_cols, "met": b.metab_cols, "mic": b.microbiome_cols},
        data=b.df,
        confounders=["age", "sex"],
        method="hima",
        n_boot=60,
        random_state=0,
    )
    assert res.method == "hima"
    assert "selected_mediators" in res.extras
    assert np.sign(res.total_effect.estimate) == np.sign(b.ground_truth_total)
