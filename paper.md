---
title: 'skinflame: Multi-omic causal mediation analysis for atopic dermatitis and inflammation research'
tags:
  - Python
  - causal inference
  - mediation analysis
  - multi-omics
  - atopic dermatitis
  - inflammation
  - bioinformatics
authors:
  - name: Vishnupriya Kannan
    orcid: 0009-0007-1326-4549
    affiliation: 1
affiliations:
  - name: Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore
    index: 1
date: 4 May 2026
bibliography: paper.bib
---

# Summary

`skinflame` is a Python library for multi-omic causal mediation analysis,
designed with atopic dermatitis (AD) and broader inflammation research in
mind. Given an exposure $X$ (e.g., a *FLG* loss-of-function variant, an
environmental pollutant, a biologic such as dupilumab), an outcome $Y$ (e.g.,
SCORAD, EASI, an itch NRS, a mood instrument), one or more high-dimensional
mediator blocks $M_1, \ldots, M_k$ (transcriptome, metabolome, microbiome,
methylome), and a set of confounders $C$, the package decomposes the total
exposure-outcome effect into direct and natural-indirect components, returns
bootstrap confidence intervals, and selects sparse sets of active mediators
with FDR control. It also accepts a user-supplied causal DAG and derives a
valid back-door (or, when needed, front-door) adjustment set, removing a
common source of bias in observational AD studies.

# Statement of Need

Mediation analysis is increasingly central to translational AD and
neuro-inflammation research, where the question is rarely *whether* a
treatment or genotype affects severity but *through which biological
pathways*. Existing software is fragmented: `mediation` [@imai2010] and
`HIMA` [@zhang2016hima] in R cover the classical and high-dimensional
cases, while `bama` and `mma` cover Bayesian and multi-mediator extensions
respectively. There is no equivalent Python tool that bundles classical
mediation, high-dimensional selection, pathway-level aggregation, and
DAG-aware adjustment behind a single API, and none of the existing tools
ship with AD-specific defaults or a working demo on AD omics data.

`skinflame` fills this gap. It targets methodologists, computational
biologists, and clinical researchers who want to ask *"how much of the
filaggrin → SCORAD effect runs through Th2 cytokines?"* in a few lines of
Python, and to do so with sound multiple-testing control and a
reproducible Streamlit demo.

# Methods

For the linear-Gaussian case the natural indirect effect is the sum of
products $\sum_j \alpha_j \beta_j$, where $\alpha_j$ is the
exposure-on-mediator coefficient and $\beta_j$ is the
mediator-on-outcome coefficient (adjusted for the exposure and other
mediators). `skinflame` implements three estimators:

1. **Classic.** OLS on each mediator regression and a joint outcome
   regression, with a paired bootstrap over the entire pipeline and BCa
   intervals.
2. **HIMA** [@zhang2016hima]. Sure-independence screening to reduce the
   mediator pool to $\sim 2n / \log n$ candidates, an $\ell_1$-penalised
   outcome model to enforce sparsity, joint-significance testing for each
   surviving mediator, and Benjamini-Hochberg FDR control over those
   joint p-values.
3. **Pathway-level.** Aggregate per-gene mediators into pathway scores
   (PCA-1 or mean-z), then run the classic estimator at the pathway
   level. Useful when biological interpretation matters more than
   resolving individual genes.

DAG-aware adjustment uses Pearl's back-door criterion [@pearl2009] over a
user-supplied `networkx` DAG, automatically excluding mediators from any
candidate adjustment set. A front-door fallback is provided when no
back-door set exists.

# Validation

Unit tests confirm that all three estimators recover known ground-truth
indirect effects on a simulated AD-flavoured cohort (sign and 95% CI
coverage). On the bundled simulator (n = 400, p = 200 transcript +
50 metabolite + 30 microbiome features, 11 active mediators), HIMA
recovered 7 of the 8 truly active mediators at FDR 0.05 and the
classical pipeline recovered the total / direct / indirect decomposition
within 5 % of ground truth.

# Software availability

- Source code: https://github.com/vp-thehuman/skinflame
- PyPI: `pip install skinflame`
- Documentation: README and `docs/` in the repository
- License: MIT

# Acknowledgements

The author thanks the open-source maintainers of `numpy`, `pandas`,
`scipy`, `scikit-learn`, `networkx`, `statsmodels`, `streamlit`, and
`GEOparse`, on which `skinflame` is built.

# References
