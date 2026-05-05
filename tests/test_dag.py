import networkx as nx
import pytest

from skinflame.dag import adjustment_set, front_door_set, validate_dag


def _classic_confounder_dag():
    """X <- C -> Y, X -> Y, X -> M -> Y."""
    g = nx.DiGraph()
    g.add_edges_from(
        [
            ("C", "X"),
            ("C", "Y"),
            ("X", "Y"),
            ("X", "M"),
            ("M", "Y"),
        ]
    )
    return g


def test_validate_dag_rejects_cycle():
    g = nx.DiGraph()
    g.add_edges_from([("A", "B"), ("B", "A")])
    with pytest.raises(ValueError):
        validate_dag(g)


def test_adjustment_set_picks_confounder_excludes_mediator():
    g = _classic_confounder_dag()
    z = adjustment_set(g, exposure="X", outcome="Y")
    assert "C" in z
    assert "M" not in z  # mediator must NOT be in adjustment set


def test_adjustment_set_must_include_respected():
    g = _classic_confounder_dag()
    g.add_node("Q")  # disconnected
    z = adjustment_set(g, exposure="X", outcome="Y", must_include=["Q"])
    assert "Q" in z


def test_front_door_when_unobserved_confounder():
    """X <- U -> Y (U unobserved), X -> M -> Y. Should yield M as front-door."""
    g = nx.DiGraph()
    g.add_edges_from([("U", "X"), ("U", "Y"), ("X", "M"), ("M", "Y")])
    fd = front_door_set(g, exposure="X", outcome="Y")
    assert fd == {"M"}
