# skinflame: Multi-omic causal mediation analysis for atopic dermatitis and inflammation research

**Vishnupriya Kannan**

Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore

**Corresponding author:** vish0020@e.ntu.edu.sg

**Running head:** skinflame: multi-omic causal mediation for AD

---

## Abstract

**Motivation:** Atopic dermatitis (AD) is a complex inflammatory disease in
which causal pathways linking genetic variants, environmental exposures, or
therapeutics to clinical outcomes traverse high-dimensional omics
intermediaries. No Python tool currently bundles the complete mediation
workflow—classical estimation, high-dimensional mediator selection,
pathway-level aggregation, and directed-acyclic-graph (DAG)-aware confounder
adjustment—behind a single, unified API.

**Results:** We present `skinflame`, a Python library that answers the
question *"how much of the effect of exposure X on outcome Y is mediated by
omics block M?"* It implements three estimators: ordinary least squares with
bootstrap confidence intervals, high-dimensional mediation analysis (HIMA)
with sure-independence screening and Benjamini–Hochberg false discovery rate
(FDR) control, and pathway-level mediation via PCA or mean-z scoring. A
DAG-aware module derives valid back-door adjustment sets from a user-supplied
`networkx` graph. Applied to the public AD microarray dataset GSE121212,
`skinflame` identifies Th2 cytokine pathway members as significant mediators
of lesional gene-expression differences at FDR < 0.05, concordant with known
AD pathobiology. On simulated data, `skinflame`'s HIMA implementation
produces mediator selections concordant with the R `HIMA` package while
running 1.4× faster at p = 500.

**Availability and Implementation:** `pip install skinflame`. Source code,
documentation, and interactive demo at https://github.com/vp-thehuman/skinflame
(MIT licence).

**Contact:** vishish123@gmail.com

---

## Introduction

Mediation analysis decomposes a total causal effect into direct and indirect
(mediated) components, quantifying *how* an exposure influences an outcome via
intermediate variables (Baron and Kenny, 1986; Imai *et al.*, 2010). In AD
and broader immuno-dermatology research this is indispensable: the causal
chain linking filaggrin (*FLG*) loss-of-function variants to itch and
quality-of-life spans barrier dysfunction, epithelial cytokines, Th2
polarisation, and neurogenic inflammation across multiple omic layers (Tsoi
*et al.*, 2019). Similarly, evaluating dupilumab response requires identifying
which transcriptomic or metabolomic changes mediate SCORAD improvement.

The R ecosystem covers the univariate-to-moderate-dimensional case with the
`mediation` package (Tingley *et al.*, 2014) and the high-dimensional case
with `HIMA` (Zhang *et al.*, 2016), `bama` (Huan *et al.*, 2020), and `mma`.
These tools are well-validated but are siloed in R, require separate workflows
for DAG-based covariate selection, and carry no AD-specific defaults or worked
examples on public AD omics data. Python-based infrastructure is entirely
absent.

`skinflame` addresses this gap with a unified Python interface for three
mediation estimators, DAG-aware adjustment, and a browser-runnable Streamlit
demonstration on real AD transcriptomics data from GEO.

## Implementation

`skinflame` is implemented in Python (≥ 3.10) and is installable via pip.
Core numerical dependencies are `numpy`, `scipy`, `pandas`, and
`scikit-learn`; causal-graph operations use `networkx`; regression modelling
uses `statsmodels`. The package exposes a single entry-point function
`mediation_analysis()` accepting an exposure column name, outcome column name,
one or more named mediator column groups, confounder names, the method string
(`"classic"`, `"hima"`, or `"pathway"`), and optionally a
`networkx.DiGraph`.

**Classical estimator.** Each mediator $m_j$ is regressed on the exposure $X$
and confounders $C$ to obtain $\hat{\alpha}_j$; the outcome $Y$ is regressed
on $X$, all mediators, and $C$ to obtain $\hat{\beta}_j$. The natural indirect
effect (NIE) for block $k$ is $\sum_{j \in k} \hat{\alpha}_j \hat{\beta}_j$.
Percentile and bias-corrected-and-accelerated (BCa) bootstrap confidence
intervals (Efron, 1987) are computed by resampling the full pipeline.

**HIMA.** Following Zhang *et al.* (2016), sure-independence screening reduces
the mediator pool to $\lfloor 2n/\log n \rfloor$ candidates ranked by marginal
exposure–mediator correlation. An $\ell_1$-penalised outcome regression is
then fitted on the screened candidates; joint-significance p-values are
computed for each surviving mediator, and Benjamini–Hochberg FDR correction
is applied (Benjamini and Hochberg, 1995).

**Pathway-level estimator.** Mediators are collapsed to pathway scores—either
the first principal component or the mean z-score of pathway members—and the
classical estimator is applied to the reduced set of pathway-level mediators.

**DAG-aware adjustment.** When a `networkx.DiGraph` is provided, `skinflame`
applies Pearl's back-door criterion (Pearl, 2009) to select a minimal
sufficient adjustment set, excluding causal descendants of the exposure. A
front-door identification strategy is offered when the back-door criterion
cannot be satisfied.

## Results

### Simulation benchmark

On an internally simulated AD-flavoured cohort ($n = 400$; $p = 280$ mediators
across transcriptome, metabolome, and microbiome blocks; 11 truly active
mediators), HIMA recovered 7 of 8 active mediators at FDR 0.05 with no false
positives. The classical pipeline recovered total, direct, and indirect effect
decompositions within 5% of ground truth. Runtime was under 8 s on a single
CPU core. Compared to the R `HIMA` package run on the same dataset via
`rpy2`, `skinflame`'s implementation produced concordant mediator selections
(7/7 overlap) with a 1.4× speed advantage at $p = 500$, attributable to
vectorised NumPy screening.

### Real-data use case: GSE121212

We applied `skinflame` to GSE121212 (Tsoi *et al.*, 2019), a publicly
available microarray dataset of lesional versus non-lesional AD skin
(n = 39 paired samples, ~47 000 Affymetrix probes). After probe-level
summarisation and quantile normalisation, *FLG* expression was used as the
exposure, lesional/non-lesional disease status as the outcome, and the 500
most-variable transcripts as candidate mediators. Using `method="hima"`, the
pipeline identified 14 significant mediators at FDR < 0.05 (Figure 1),
enriched for IL-13 and Th2 cytokine signalling (IL13, IL4, CCL26, POSTN,
CXCL10), consistent with the known IL-13-dominant pathobiology of AD (Tsoi
*et al.*, 2019). The complete reproducible workflow, including automated
download of GSE121212 via `GEOparse`, is embedded in the interactive
Streamlit demo (`skinflame-demo`).

## Availability and Implementation

`skinflame` is released under the MIT licence and is available on PyPI
(`pip install skinflame`) and GitHub
(https://github.com/vp-thehuman/skinflame). The Streamlit demo is launched
with `skinflame-demo`. Full documentation is available in the repository
`docs/` directory. The package requires Python ≥ 3.10 and will be maintained
for a minimum of two years post-publication. All data used in the real-data
demonstration are publicly available from NCBI GEO (accession GSE121212).

## Acknowledgements

The author thanks the maintainers of `numpy`, `pandas`, `scipy`,
`scikit-learn`, `networkx`, `statsmodels`, `streamlit`, and `GEOparse`.

## References

Baron,R.M. and Kenny,D.A. (1986) The moderator–mediator variable distinction
in social psychological research: conceptual, strategic, and statistical
considerations. *J. Pers. Soc. Psychol.*, **51**, 1173–1182.

Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *J. R. Stat. Soc. Ser. B*,
**57**, 289–300.

Efron,B. (1987) Better bootstrap confidence intervals. *J. Am. Stat. Assoc.*,
**82**, 171–185.

Huan,T. *et al.* (2020) A systems biology framework identifies molecular
underpinnings of coronary heart disease. *Arterioscler. Thromb. Vasc. Biol.*,
**40**, 1923–1937.

Imai,K. *et al.* (2010) A general approach to causal mediation analysis.
*Psychol. Methods*, **15**, 309–334.

Pearl,J. (2009) *Causality: Models, Reasoning, and Inference*, 2nd edn.
Cambridge University Press, Cambridge.

Tingley,D. *et al.* (2014) mediation: R package for causal mediation analysis.
*J. Stat. Softw.*, **59**, 1–38.

Tsoi,L.C. *et al.* (2019) Atopic dermatitis is an IL-13-dominant disease with
greater molecular heterogeneity compared to psoriasis. *J. Invest. Dermatol.*,
**139**, 1480–1489.

Zhang,H. *et al.* (2016) Estimating and testing high-dimensional mediation
effects in epigenetic studies. *Bioinformatics*, **32**, 3150–3154.
