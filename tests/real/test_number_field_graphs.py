"""Real graph comparison respects the chosen embedding of the constant field."""

import pytest
import sympy as sp

from alggeom.decomposition.rational import RationalDecomposition
from alggeom.decomposition.real import verify_decomposition

x, t, y = sp.symbols("x t y", real=True)


@pytest.mark.parametrize("alpha", [sp.sqrt(2), -sp.sqrt(2), sp.real_root(2, 3)])
def test_real_embeddings(alpha):
    constants = sp.QQ.algebraic_field(alpha)
    result = RationalDecomposition(x**4 + alpha, (x,), t, constants, t**2 + alpha, x**2)
    comparison = verify_decomposition(result, output=y)
    assert comparison.equivalent and comparison.certified and comparison.function_verified


def test_nonreal_embedding():
    constants = sp.QQ.algebraic_field(sp.sqrt(-3))
    result = RationalDecomposition(x**4 - 3, (x,), t, constants, t**2 - 3, x**2)
    with pytest.raises(ValueError, match="real embedding"):
        verify_decomposition(result, output=y)
