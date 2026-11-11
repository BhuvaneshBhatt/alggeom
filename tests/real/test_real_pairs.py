"""Real domains, poles and branch choices for two-component certificates."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.decomposition.algebraic import AlgebraicDecomposition
from alggeom.decomposition.discovery import search_algebraic
from alggeom.decomposition.rational import RationalDecomposition
from alggeom.decomposition.real import verify_decomposition

x, y, z = sp.symbols("x y z", real=True)


def test_automatic_branch():
    search = search_algebraic(x * x - z * z - 1, (x, z), parameter=y)
    assert search.complete and search.verify()
    (D,) = search.decompositions
    result = verify_decomposition(D, target_selector=x >= 0, outer_selector=x >= 0)
    assert result.equivalent and result.certified and result.function_verified
    wrong = verify_decomposition(D, target_selector=x >= 0, outer_selector=x <= 0)
    assert not wrong.equivalent


def test_rational_graph():
    D = RationalDecomposition(x**4, (x,), y, sp.QQ, y * y, x * x)
    R = verify_decomposition(D, output=z)
    assert R.equivalent and R.certified and R.function_verified
    assert not verify_decomposition(D, output=z, inner_selector=y < 0).equivalent


def test_internal_pole():
    D = RationalDecomposition(1 / (x**4 + 1), (x,), y, sp.QQ, y * y / (y * y + 1), 1 / x**2)
    R = verify_decomposition(D, output=z)
    assert not R.equivalent and R.certified and (not R.function_verified)
    restricted = verify_decomposition(D, output=z, domain=sp.Ne(x, 0))
    assert restricted.equivalent and restricted.certified and restricted.function_verified


def test_original_hole():
    original = sp.Mul(x**4 * (x - 1), sp.Pow(x - 1, -1, evaluate=False), evaluate=False)
    D = RationalDecomposition(original, (x,), y, sp.QQ, y * y, x * x)
    assert D.verify()
    assert not verify_decomposition(D, output=z).equivalent
    assert verify_decomposition(D, output=z, domain=sp.Ne(x, 1)).function_verified


def test_principal_root():
    D = AlgebraicDecomposition(
        sp.Poly(x**4 - z, x, z, domain=sp.QQ),
        sp.Poly(x * x - y, x, y, domain=sp.QQ),
        sp.Poly(y * y - z, y, z, domain=sp.QQ),
        1,
        sp.Integer(1),
    )
    assert D.verify()
    relation = verify_decomposition(D, domain=z >= 0)
    assert relation.equivalent and relation.certified and (not relation.function_verified)
    principal = verify_decomposition(
        D, domain=z >= 0, target_selector=x >= 0, outer_selector=x >= 0, inner_selector=y >= 0
    )
    assert principal.equivalent and principal.certified and principal.function_verified
    wrong = verify_decomposition(
        D, domain=z > 0, target_selector=x >= 0, outer_selector=x >= 0, inner_selector=y <= 0
    )
    assert not wrong.equivalent


def test_invalid_certificate():
    D = RationalDecomposition(x**4, (x,), y, sp.QQ, y * y, x * x)
    with pytest.raises(ValueError, match="failed verification"):
        verify_decomposition(replace(D, outer=y * y + 1), output=z)
    with pytest.raises(ValueError):
        verify_decomposition(D, output=x)
    with pytest.raises(ValueError, match="undeclared"):
        verify_decomposition(D, output=z, inner_selector=z > 0)
    with pytest.raises(ValueError, match="certificate"):
        verify_decomposition(object())


def test_complex_constants():
    K = sp.QQ.algebraic_field(sp.I)
    D = RationalDecomposition(x**4 + sp.I, (x,), y, K, y * y + sp.I, x * x)
    with pytest.raises(ValueError, match="real embedding"):
        verify_decomposition(D, output=z)
