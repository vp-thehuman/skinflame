# skinflame

[![PyPI](https://img.shields.io/pypi/v/skinflame.svg)](https://pypi.org/project/skinflame/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20034163.svg)](https://doi.org/10.5281/zenodo.20034163)

**Multi-omic causal mediation analysis for atopic dermatitis (AD) and inflammation research.**

`skinflame` answers questions of the form:

> *How much of the effect of exposure E (e.g. a filaggrin variant, a pollutant,
> dupilumab) on outcome Y (SCORAD, EASI, itch NRS, mood) flows through omics
> block M (transcriptome, metabolome, microbiome, methylome)?*

It is designed for the **AD ↔ mood / inflammation** research interface, where
mediator blocks are high-dimensional and frequently nested.

## Features

- **Classical mediation** — total / direct / natural-indirect effects with
  product-of-coefficients, Sobel, joint-significance, and percentile / BCa
  bootstrap confidence intervals.
- **HIMA — high-dimensional mediation** — sure-independence screening
  followed by MCP-style penalisation on the joint significance of the
  exposure→mediator→outcome paths, with FDR control.
- **Pathway-level mediation** — aggregate gene/metabolite mediators into
  pathway scores (PCA-1 or mean-z) and test the natural indirect effect at the
  pathway level.
- **DAG-aware adjustment** — pass a `networkx.DiGraph` and `skinflame` derives
  a valid back-door adjustment set (or warns if only front-door is available).
- **Streamlit demo** — explore the methods on a synthetic AD-flavoured cohort
  or pull GSE121212 (real AD lesional vs non-lesional transcriptomics) live
  from GEO.

## Install

```bash
pip install skinflame                # core
pip install "skinflame[demo]"        # + Streamlit + plotly
pip install "skinflame[geo]"         # + GEOparse for GSE121212
pip install "skinflame[all]"         # everything
```

## Quickstart

```python
import numpy as np
from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort

data = simulate_ad_cohort(n=400, p_transcript=200, p_metab=50, seed=0)

result = mediation_analysis(
    exposure="filaggrin_LoF",
    outcome="SCORAD",
    mediators={"transcriptome": data.transcript_cols,
               "metabolome":   data.metab_cols},
    confounders=["age", "sex", "BMI"],
    data=data.df,
    n_boot=1000,
    method="hima",          # "classic" | "hima" | "pathway"
    random_state=0,
)

print(result.summary())
```

Outputs include `total_effect`, `direct_effect`, `indirect_effect`, the
proportion mediated, per-block decomposition, and (with HIMA) a
`selected_mediators` table with FDR-adjusted joint p-values.

## Streamlit demo

```bash
skinflame-demo            # or: streamlit run demo/streamlit_app.py
```

Toggle between synthetic data (default, always works) and GSE121212
(real microarray, downloaded on first run via `GEOparse`).

## Why this exists

Multi-omic mediation is the empirical bridge between *what we measure*
(SNPs, exposures, drugs) and *what patients feel* (itch, sleep, depression).
Existing tooling is fragmented across R packages (`mediation`, `HIMA`,
`bama`, `mma`) with little Python coverage and almost no AD-specific
defaults. `skinflame` aims to be the obvious Python entry point.

## Citation

If you use `skinflame` in academic work, please cite:

> Vishnupriya, K. (2026). *skinflame: Multi-omic causal mediation analysis for
> atopic dermatitis and inflammation research.* [submission pending].

## License

MIT — see [LICENSE](LICENSE).
