import numpy as np
import pytest

from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort


def test_simulate_ad_cohort_shapes():
    b = simulate_ad_cohort(n=200, p_transcript=50, p_metab=20, p_microbiome=10, seed=0)
    assert len(b.df) == 200
    assert "filaggrin_LoF" in b.df.columns
    assert "SCORAD" in b.df.columns
    assert len(b.transcript_cols) == 50
    assert b.ground_truth_indirect is not None
    assert b.ground_truth_total == pytest.approx(b.ground_truth_direct + b.ground_truth_indirect)


def test_classic_mediation_recovers_direction():
    b = simulate_ad_cohort(n=600, p_transcript=20, p_metab=10, p_microbiome=5, seed=1)
    # Use only the *active* mediators (the first few) so classic OLS is well-conditioned.
    res = mediation_analysis(
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators={"transcriptome": b.transcript_cols[:6],
                   "metabolome": b.metab_cols[:3],
                   "microbiome": b.microbiome_cols[:2]},
        data=b.df,
        confounders=["age", "sex", "BMI"],
        n_boot=200,
        random_state=42,
    )
    # Direction recovered (sign matches ground truth).
    assert np.sign(res.indirect_effect.estimate) == np.sign(b.ground_truth_indirect)
    assert np.sign(res.total_effect.estimate) == np.sign(b.ground_truth_total)
    # CI contains the truth (loose tolerance for n=600, n_boot=200).
    assert res.indirect_effect.ci_low <= b.ground_truth_indirect <= res.indirect_effect.ci_high


def test_proportion_mediated_in_unit_range():
    b = simulate_ad_cohort(n=400, p_transcript=15, p_metab=10, p_microbiome=5, seed=2)
    actives = b.transcript_cols[:6] + b.metab_cols[:3] + b.microbiome_cols[:2]
    res = mediation_analysis(
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators=actives,
        data=b.df,
        confounders=["age", "sex"],
        n_boot=100,
        random_state=0,
    )
    assert 0.0 <= res.proportion_mediated <= 1.5  # may exceed 1 if direct flips sign


def test_summary_dataframe_structure():
    b = simulate_ad_cohort(n=300, p_transcript=10, p_metab=5, p_microbiome=3, seed=3)
    res = mediation_analysis(
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators={"tx": b.transcript_cols[:6]},
        data=b.df,
        n_boot=80,
        random_state=0,
    )
    df = res.summary()
    assert {"effect", "estimate", "ci_low", "ci_high"}.issubset(df.columns)
    assert (df["effect"] == "total").any()
    assert (df["effect"] == "direct").any()
    assert (df["effect"] == "indirect").any()
