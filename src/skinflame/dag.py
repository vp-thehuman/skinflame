"""DAG-aware adjustment for mediation analysis.

Given a user-supplied causal DAG (a `networkx.DiGraph`) over the variables in
the analysis, return a *valid* back-door adjustment set for the
exposure -> outcome contrast. Mediators on the X -> Y path are explicitly
excluded from the adjustment set (adjusting for a mediator blocks part of the
indirect effect we want to measure).

The implementation follows Pearl's back-door criterion: a set Z of nodes
satisfies the back-door criterion relative to (X, Y) iff
  (i) no node in Z is a descendant of X, and
  (ii) Z d-separates X and Y in the graph obtained by deleting all outgoing
       edges from X.

We additionally screen out mediators (descendants of X that are ancestors of
Y) from candidate adjustment sets, since adjusting for a mediator changes
the estimand from total to direct effect.
"""
from __future__ import annotations

from collections.abc import Iterable

import networkx as nx


def validate_dag(g: nx.DiGraph, exposure: str | None = None, outcome: str | None = None) -> None:
    """Raise if `g` is not a DAG or if exposure/outcome are missing."""
    if not isinstance(g, nx.DiGraph):
        raise TypeError("`dag` must be a networkx.DiGraph")
    if not nx.is_directed_acyclic_graph(g):
        raise ValueError("Graph contains a cycle; not a DAG.")
    for node in (exposure, outcome):
        if node is not None and node not in g:
            raise KeyError(f"Node {node!r} not in DAG.")


def _moralize(g: nx.DiGraph) -> nx.Graph:
    """Build the moral graph (used in the m-separation / back-door check)."""
    m = nx.Graph()
    m.add_nodes_from(g.nodes)
    for n in g.nodes:
        parents = list(g.predecessors(n))
        # connect each pair of parents
        for i in range(len(parents)):
            for j in range(i + 1, len(parents)):
                m.add_edge(parents[i], parents[j])
        for p in parents:
            m.add_edge(p, n)
    return m


def d_separated(
    g: nx.DiGraph, x: Iterable[str] | str, y: Iterable[str] | str, z: Iterable[str] | None = None
) -> bool:
    """Pearl-style d-separation between sets X and Y given Z.

    Uses networkx's built-in implementation when available (>=3.1), otherwise
    falls back to the moralization-based ancestral-graph check.
    """
    X = {x} if isinstance(x, str) else set(x)
    Y = {y} if isinstance(y, str) else set(y)
    Z = set(z) if z else set()
    try:
        # networkx >= 3.5 uses is_d_separator; older versions use d_separated
        if hasattr(nx, "is_d_separator"):
            return nx.is_d_separator(g, X, Y, Z)
        return nx.algorithms.d_separation.d_separated(g, X, Y, Z)
    except Exception:
        # Fallback: ancestral moral graph check
        ancestors = set().union(*(nx.ancestors(g, n) | {n} for n in X | Y | Z))
        sub = g.subgraph(ancestors).copy()
        moral = _moralize(sub)
        moral.remove_nodes_from(Z - X - Y)
        if not (X & set(moral.nodes)) or not (Y & set(moral.nodes)):
            return True
        return not any(
            nx.has_path(moral, a, b) for a in X if a in moral for b in Y if b in moral
        )


def _candidates(
    g: nx.DiGraph,
    exposure: str,
    outcome: str,
    forbidden: set[str],
) -> set[str]:
    """Variables eligible to be in an adjustment set:
    not descendants of exposure, not exposure/outcome themselves, not forbidden.
    """
    desc = nx.descendants(g, exposure) | {exposure}
    return set(g.nodes) - desc - {outcome} - set(forbidden)


def adjustment_set(
    g: nx.DiGraph,
    exposure: str,
    outcome: str,
    must_include: Iterable[str] | None = None,
    forbidden: Iterable[str] | None = None,
) -> set[str]:
    """Return a valid back-door adjustment set Z for (exposure -> outcome).

    Strategy: start with all back-door-eligible candidates (non-descendants
    of exposure), verify the back-door criterion via d-separation in the
    graph with all outgoing edges from `exposure` removed; if that succeeds,
    greedily prune nodes that aren't needed.

    Falls back to ``set(must_include)`` if no valid set is found (and warns).

    Parameters
    ----------
    must_include : variables to force into the set.
    forbidden    : variables to exclude (e.g. known mediators, post-treatment
                   variables). Mediators on X→Y are auto-detected and added
                   to ``forbidden``.
    """
    validate_dag(g, exposure=exposure, outcome=outcome)
    must_include = set(must_include or [])
    forbidden = set(forbidden or [])

    # Auto-detect mediators: descendants of X that are ancestors of Y
    desc_x = nx.descendants(g, exposure)
    anc_y = nx.ancestors(g, outcome)
    mediators = desc_x & anc_y
    forbidden |= mediators

    # Build mutilated graph (remove edges out of exposure) for back-door check
    g_bd = g.copy()
    g_bd.remove_edges_from(list(g.out_edges(exposure)))

    cands = _candidates(g, exposure, outcome, forbidden)

    # Try the maximal eligible set first
    Z = set(cands) | must_include
    if not d_separated(g_bd, exposure, outcome, Z):
        # No valid set exists with this candidate pool
        import warnings

        warnings.warn(
            "No valid back-door adjustment set found; returning `must_include` only. "
            "Consider front-door adjustment or revisiting the DAG.",
            stacklevel=2,
        )
        return set(must_include)

    # Greedy prune: drop any node we can lose while preserving d-separation
    pruned = set(Z)
    # Sort for determinism; prefer keeping `must_include`
    for node in sorted(Z - must_include):
        trial = pruned - {node}
        if d_separated(g_bd, exposure, outcome, trial):
            pruned = trial
    return pruned


def front_door_set(g: nx.DiGraph, exposure: str, outcome: str) -> set[str] | None:
    """Return a set M satisfying Pearl's front-door criterion, or None.

    Conditions:
      1. M intercepts every directed path from exposure to outcome.
      2. There is no unblocked back-door path from exposure to M.
      3. Every back-door path from M to outcome is blocked by exposure.
    """
    validate_dag(g, exposure=exposure, outcome=outcome)

    # Candidate M: descendants of X that are ancestors of Y
    desc_x = nx.descendants(g, exposure)
    anc_y = nx.ancestors(g, outcome)
    cand = sorted(desc_x & anc_y)

    for k in range(1, len(cand) + 1):
        from itertools import combinations

        for subset in combinations(cand, k):
            M = set(subset)
            # (1) Every directed X→Y path passes through M: remove M, no path remains.
            sub = g.copy()
            sub.remove_nodes_from(M)
            if nx.has_path(sub, exposure, outcome):
                continue
            # (2) No back-door path X→M unblocked (empty conditioning set).
            g_bd = g.copy()
            g_bd.remove_edges_from(list(g.out_edges(exposure)))
            if not d_separated(g_bd, exposure, M, set()):
                continue
            # (3) Every back-door path M→Y blocked by X.
            g_bd_m = g.copy()
            for m in M:
                g_bd_m.remove_edges_from(list(g.out_edges(m)))
            if not d_separated(g_bd_m, M, outcome, {exposure}):
                continue
            return M
    return None
