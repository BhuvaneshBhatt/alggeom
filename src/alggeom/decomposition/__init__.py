"""Exact compositions of rational maps and implicit algebraic functions."""

from ._api import (
    common_right_component,
    decompose,
    decompose_relation,
    decomposition_condition,
    decompositions,
    real_equivalent,
)
from ._validation import DecompositionError, DegreeLimitError, IncompleteDecompositionError

__all__ = [
    "common_right_component",
    "decompose",
    "decompositions",
    "decompose_relation",
    "decomposition_condition",
    "real_equivalent",
    "DecompositionError",
    "DegreeLimitError",
    "IncompleteDecompositionError",
]
