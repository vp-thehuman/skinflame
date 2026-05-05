
from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort


def test_pathway_first_pathway_is_significant():
    b = simulate_ad_cohort(n=400, p_transcript=80, p_metab=10, p_microbiome=5, seed=0)
    res = mediation_analysis(
        exposure="filaggrin_LoF",
        outcome="SCORAD",
        mediators=b.transcript_cols,  # ignored when method=pathway
        data=b.df,
        method="pathway",
        pathways=b.pathway_map,
        confounders=["age", "sex"],
        n_boot=60,
        random_state=0,
    )
    tbl = res.extras["pathway_table"]
    assert "PATHWAY_01" in tbl["pathway"].tolist()
    # The pathway containing the active mediators should rank #1 and be FDR-significant.
    assert tbl.iloc[0]["pathway"] == "PATHWAY_01"
    assert bool(tbl.iloc[0]["fdr_significant"])
