"""skinflame — multi-omic causal mediation analysis for AD / inflammation research."""

from .core import (
    MediationResult,
    bootstrap_mediation,
    mediation_analysis,
)
from .dag import adjustment_set, validate_dag
from .hima import hima
from .pathway import pathway_mediation

__version__ = "0.1.0"

__all__ = [
    "MediationResult",
    "mediation_analysis",
    "bootstrap_mediation",
    "hima",
    "pathway_mediation",
    "adjustment_set",
    "validate_dag",
    "__version__",
]
