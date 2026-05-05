# Theory & methods

## 1. The mediation estimand

For a binary or continuous exposure $X$, an outcome $Y$, a (possibly
multivariate) mediator $M$, and confounders $C$, the natural direct effect
(NDE) and natural indirect effect (NIE) are defined in the potential-outcomes
framework (Imai, Keele & Tingley, 2010):

$$
\mathrm{NDE}(x, x^\*) = \mathbb{E}\!\left[Y(x, M(x^\*)) - Y(x^\*, M(x^\*)) \mid C\right],
$$

$$
\mathrm{NIE}(x, x^\*) = \mathbb{E}\!\left[Y(x, M(x)) - Y(x, M(x^\*)) \mid C\right].
$$

Their sum is the total effect (TE). For the linear-Gaussian model

$$
M_j = \alpha_{0j} + \alpha_j X + \boldsymbol{\gamma}_M^\top C + e_{Mj},\quad
Y = \beta_0 + \beta_X X + \boldsymbol{\beta}_M^\top M + \boldsymbol{\gamma}_Y^\top C + e_Y,
$$

these collapse to the familiar Baron-Kenny decomposition,

$$
\mathrm{NIE} = \sum_j \alpha_j \beta_{M,j},\qquad \mathrm{NDE} = \beta_X.
$$

`skinflame.core` estimates these by OLS and bootstraps the entire pipeline
(both the mediator regressions and the outcome regression) jointly, which
yields valid CIs even when $\sum_j \alpha_j \beta_{M,j}$ is the sum of many
small terms.

## 2. High-dimensional mediation (HIMA)

When the mediator block is genome-scale, classical OLS is rank-deficient. We
follow Zhang et al. (2016):

1. **Sure-Independence Screening (SIS).** Rank mediators by partial
   correlation $|\mathrm{corr}(M_j, Y \mid X, C)|$ and keep the top
   $d \approx 2n / \log n$.
2. **Penalized regression.** Fit $Y \sim X + M_{\text{screened}} + C$ with a
   sparsity-inducing penalty (we use `LassoCV` as an open-source
   stand-in for MCP/SCAD).
3. **Joint significance.** For each surviving mediator, compute
   $p_\alpha$ from $X \to M_j$ and $p_\beta$ from $M_j \to Y$ (in the
   multi-mediator outcome model). The joint p-value is
   $p_j^{JS} = \max(p_\alpha, p_\beta)$ — this is the maximum of two
   correlated tests and is conservative.
4. **FDR control.** Benjamini-Hochberg on $\{p_j^{JS}\}$ at level $q$
   (default 0.05).

## 3. Pathway-level mediation

Per-gene mediation suffers from low power and replication issues. We
aggregate genes into pathway scores via PCA-1 (sign-aligned to the column
mean) or simple mean-z, then apply the classic mediation pipeline at the
pathway level. This trades a small amount of resolution for a large gain
in interpretability and statistical power, and gives results that map
directly onto Reactome / KEGG / Hallmark gene sets.

## 4. DAG-aware adjustment

`skinflame.dag.adjustment_set` accepts a user-supplied causal DAG and
returns a *valid back-door adjustment set* — that is, a set $Z$ such that:

  - no element of $Z$ is a descendant of $X$ (no post-treatment bias), and
  - $Z$ d-separates $X$ from $Y$ in the graph with all out-edges of $X$
    removed.

Mediators on $X \to Y$ are auto-detected and excluded (adjusting for them
would change the estimand from total to direct). When no valid back-door
set exists, the function emits a warning and the user can fall back to
`front_door_set`, which implements Pearl's front-door criterion.

## 5. Bootstrap details

Both percentile and BCa bootstrap CIs are supported. BCa uses a
jackknife-on-bootstrap-samples approximation for the acceleration term
(a common pragmatic shortcut when refitting the original statistic per
jackknife replicate is too expensive). For HIMA, the bootstrap is applied
*after* selection — so reported intervals reflect the post-selection
estimator and should not be interpreted as standard frequentist CIs for the
underlying $\alpha_j \beta_j$ products. Joint-significance + FDR is the
inferential workhorse for individual mediators.
