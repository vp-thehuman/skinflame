"""Streamlit demo for `skinflame`.

Run with:
    streamlit run demo/streamlit_app.py
or, after pip-install:
    skinflame-demo
"""
from __future__ import annotations

import io
import textwrap

import networkx as nx
import numpy as np
import pandas as pd
import streamlit as st

from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort
from skinflame.plots import effects_dataframe

st.set_page_config(page_title="skinflame — multi-omic mediation", layout="wide")

st.title("skinflame")
st.caption(
    "Multi-omic causal mediation analysis for atopic dermatitis & inflammation research."
)

# ---------------------------------------------------------------------------
# Sidebar: dataset choice + analysis controls
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("1. Data")
    source = st.radio(
        "Dataset",
        ["Synthetic AD cohort (default)", "GSE121212 (real microarray)"],
        index=0,
        help=(
            "Synthetic data has known ground truth for mediation effects. "
            "GSE121212 requires `pip install skinflame[geo]` and downloads "
            "~50 MB on first use."
        ),
    )
    if source.startswith("Synthetic"):
        n = st.slider("Sample size", 100, 1500, 400, step=50)
        p_tx = st.slider("# transcriptome mediators", 20, 500, 200, step=20)
        p_met = st.slider("# metabolome mediators", 5, 200, 50, step=5)
        seed = st.number_input("Random seed", value=0, step=1)
    else:
        st.info(
            "GSE121212 = AD lesional vs non-lesional skin (Tsoi et al., 2019). "
            "First load is slow (downloads from GEO)."
        )
        top_genes = st.slider("Top-variance genes", 100, 2000, 500, step=100)

    st.header("2. Method")
    method = st.selectbox(
        "Mediation method",
        ["classic", "hima", "pathway"],
        index=1,
        help=(
            "classic: OLS + bootstrap. "
            "hima: SIS + Lasso + joint-significance + FDR for high-dim mediators. "
            "pathway: PCA-aggregate per pathway, then classic."
        ),
    )
    n_boot = st.slider("Bootstrap replicates", 100, 5000, 500, step=100)
    ci_kind = st.selectbox("CI type", ["percentile", "bca"], index=0)
    use_dag = st.checkbox(
        "Use built-in DAG for adjustment",
        value=True,
        help="The synthetic generator's true DAG is provided; toggling derives a back-door set.",
    )
    run = st.button("Run analysis", type="primary", use_container_width=True)


# ---------------------------------------------------------------------------
# Build dataset on demand
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _make_synth(n, p_tx, p_met, seed):
    return simulate_ad_cohort(
        n=int(n),
        p_transcript=int(p_tx),
        p_metab=int(p_met),
        p_microbiome=30,
        seed=int(seed),
    )


@st.cache_data(show_spinner=True)
def _load_geo(top_genes):
    from skinflame.data import load_gse121212

    return load_gse121212(top_var_genes=int(top_genes))


if source.startswith("Synthetic"):
    bundle = _make_synth(n, p_tx, p_met, seed)
    exposure = "filaggrin_LoF"
    outcome = "SCORAD"
    mediators = {
        "transcriptome": bundle.transcript_cols,
        "metabolome": bundle.metab_cols,
        "microbiome": bundle.microbiome_cols,
    }
    confounders = bundle.confounder_cols
    pathways = bundle.pathway_map
else:
    try:
        bundle = _load_geo(top_genes)
        exposure = "lesional"
        outcome = "severity_score"
        mediators = {"transcriptome": bundle.transcript_cols}
        confounders = []
        pathways = None
    except Exception as e:  # pragma: no cover
        st.error(f"Could not load GSE121212: {e}")
        st.stop()

c1, c2 = st.columns([1, 1])
with c1:
    st.subheader("Sample")
    st.write(f"n = {len(bundle.df)}, p_mediators = {sum(len(v) for v in mediators.values())}")
    st.dataframe(bundle.df.head(10), use_container_width=True, height=240)

with c2:
    st.subheader("Schema")
    schema_md = textwrap.dedent(
        f"""
        - **Exposure**: `{exposure}`
        - **Outcome**: `{outcome}`
        - **Mediator blocks**:
        """
    )
    for blk, cols in mediators.items():
        schema_md += f"  - `{blk}` (n={len(cols)})\n"
    if confounders:
        schema_md += f"- **Confounders**: {', '.join(confounders)}\n"
    if bundle.ground_truth_total is not None:
        schema_md += "\n**Ground truth (synthetic):**\n"
        schema_md += f"  - total = {bundle.ground_truth_total:.3f}\n"
        schema_md += f"  - direct = {bundle.ground_truth_direct:.3f}\n"
        schema_md += f"  - indirect = {bundle.ground_truth_indirect:.3f}\n"
        for k, v in bundle.ground_truth_per_block.items():
            schema_md += f"    - indirect[{k}] = {v:.3f}\n"
    st.markdown(schema_md)


# ---------------------------------------------------------------------------
# DAG construction (for the synthetic case we know the true graph)
# ---------------------------------------------------------------------------
def build_synth_dag(bundle):
    g = nx.DiGraph()
    g.add_node("filaggrin_LoF")
    g.add_node("SCORAD")
    for c in bundle.confounder_cols:
        g.add_node(c)
        g.add_edge(c, "filaggrin_LoF")
        g.add_edge(c, "SCORAD")
    # Active mediators only — keeps the DAG manageable
    actives = (
        bundle.transcript_cols[:6]
        + bundle.metab_cols[:3]
        + bundle.microbiome_cols[:2]
    )
    for m in actives:
        g.add_node(m)
        g.add_edge("filaggrin_LoF", m)
        g.add_edge(m, "SCORAD")
    g.add_edge("filaggrin_LoF", "SCORAD")
    return g


dag = build_synth_dag(bundle) if (use_dag and source.startswith("Synthetic")) else None

# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------
if run:
    with st.spinner(f"Running {method} mediation with {n_boot} bootstrap replicates..."):
        kwargs = {}
        if method == "pathway":
            kwargs["pathways"] = pathways or {}
        result = mediation_analysis(
            exposure=exposure,
            outcome=outcome,
            mediators=mediators if method != "pathway" else mediators,
            data=bundle.df,
            confounders=confounders,
            method=method,
            n_boot=int(n_boot),
            ci=ci_kind,
            dag=dag,
            random_state=0,
            **kwargs,
        )

    st.success("Done.")
    st.subheader("Effects")

    summ = result.summary()
    st.dataframe(summ, use_container_width=True)

    # Forest-style plot via plotly
    try:
        import plotly.graph_objects as go

        df = effects_dataframe(result)
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=df["estimate"],
                y=df["effect"],
                mode="markers",
                error_x=dict(
                    type="data",
                    symmetric=False,
                    array=df["ci_high"] - df["estimate"],
                    arrayminus=df["estimate"] - df["ci_low"],
                ),
                marker=dict(size=10, color="#1b4965"),
            )
        )
        fig.add_vline(x=0, line_dash="dash", line_color="grey")
        fig.update_layout(
            xaxis_title="Estimate (95% CI)", yaxis_autorange="reversed", height=380
        )
        st.plotly_chart(fig, use_container_width=True)
    except Exception:
        pass

    cols = st.columns(3)
    cols[0].metric("Total effect", f"{result.total_effect.estimate:.3f}")
    cols[1].metric("Direct effect", f"{result.direct_effect.estimate:.3f}")
    cols[2].metric(
        "Proportion mediated",
        f"{(result.proportion_mediated or 0):.1%}" if result.proportion_mediated else "—",
    )

    if "selected_mediators" in result.extras:
        st.subheader("HIMA — selected mediators (FDR-controlled)")
        st.dataframe(result.extras["selected_mediators"], use_container_width=True)

    if "pathway_table" in result.extras:
        st.subheader("Pathway-level results (FDR-controlled)")
        st.dataframe(result.extras["pathway_table"], use_container_width=True)

    st.subheader("Per-mediator joint significance (top 25)")
    if result.per_mediator is not None:
        st.dataframe(
            result.per_mediator.sort_values("joint_sig_p").head(25),
            use_container_width=True,
        )

    # Download CSV
    csv = io.StringIO()
    summ.to_csv(csv, index=False)
    st.download_button(
        "Download results CSV",
        csv.getvalue(),
        file_name="skinflame_results.csv",
        mime="text/csv",
    )
else:
    st.info("Configure the run in the sidebar and click **Run analysis**.")
