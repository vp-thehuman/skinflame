# Quickstart

## Install

```bash
pip install -e ".[demo,dev]"
```

## A first analysis

```python
from skinflame import mediation_analysis
from skinflame.data import simulate_ad_cohort

bundle = simulate_ad_cohort(n=400, p_transcript=200, p_metab=50, seed=0)

result = mediation_analysis(
    exposure="filaggrin_LoF",
    outcome="SCORAD",
    mediators={
        "transcriptome": bundle.transcript_cols,
        "metabolome":   bundle.metab_cols,
        "microbiome":   bundle.microbiome_cols,
    },
    confounders=["age", "sex", "BMI"],
    data=bundle.df,
    method="hima",
    n_boot=500,
    random_state=0,
)
print(result.summary())
print(result.extras["selected_mediators"].head())
```

## With a DAG

```python
import networkx as nx

g = nx.DiGraph()
g.add_edges_from([
    ("age",       "filaggrin_LoF"),
    ("age",       "SCORAD"),
    ("sex",       "SCORAD"),
    ("BMI",       "SCORAD"),
    ("filaggrin_LoF", "T0001"),
    ("filaggrin_LoF", "T0002"),
    ("T0001",         "SCORAD"),
    ("T0002",         "SCORAD"),
    ("filaggrin_LoF", "SCORAD"),
])

result = mediation_analysis(
    exposure="filaggrin_LoF",
    outcome="SCORAD",
    mediators=["T0001", "T0002"],
    data=bundle.df,
    dag=g,
    n_boot=500,
)
```

`skinflame` derives a back-door adjustment set from the DAG (and unions it
with any `confounders=` you pass).

## Pathway-level

```python
result = mediation_analysis(
    exposure="filaggrin_LoF",
    outcome="SCORAD",
    mediators=bundle.transcript_cols,        # ignored when method="pathway"
    pathways=bundle.pathway_map,
    data=bundle.df,
    method="pathway",
    n_boot=500,
)
print(result.extras["pathway_table"])
```
