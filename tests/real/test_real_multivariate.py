"""Real graph verification with all multivariate input coordinates free."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.decomposition.rational import RationalDecomposition, search_rational
from alggeom.decomposition.real import verify_decomposition

x, y, t, z = sp.symbols("x y t z", real=True)


def test_free_coordinates():
    search = search_rational((x * x + y * y) ** 2, (x, y), parameter=t)
    result = verify_decomposition(search.decompositions[0], output=z)
    assert result.equivalent and result.certified and result.function_verified
    assert result.comparison.variables == (x, y, z)
    assert result.target_formula.free_symbols == {x, y, z}


def test_incompatible_selector():
    decomposition = RationalDecomposition(
        (x * x + y * y) ** 2, (x, y), t, sp.QQ, t * t, x * x + y * y
    )
    result = verify_decomposition(decomposition, output=z, inner_selector=t < 0)
    assert result.composition_formula is sp.false
    assert not result.equivalent and result.certified and (not result.function_verified)
    assert result.comparison.failed_direction == "lhs_implies_rhs"
    assert result.comparison.method == "quantified_symmetric_difference"


def test_original_hole():
    original = sp.Mul(
        (x * x + y * y) ** 2 * (x - 1), sp.Pow(x - 1, -1, evaluate=False), evaluate=False
    )
    decomposition = RationalDecomposition(original, (x, y), t, sp.QQ, t * t, x * x + y * y)
    unrestricted = verify_decomposition(decomposition, output=z)
    assert not unrestricted.equivalent and unrestricted.certified
    restricted = verify_decomposition(decomposition, output=z, domain=sp.Ne(x, 1))
    assert restricted.equivalent and restricted.certified and restricted.function_verified


def test_inner_poles():
    decomposition = RationalDecomposition(
        1 / (x**4 + 1), (x, y), t, sp.QQ, t * t / (t * t + 1), 1 / x**2
    )
    assert decomposition.verify()
    unrestricted = verify_decomposition(decomposition, output=z)
    assert not unrestricted.equivalent and unrestricted.certified
    restricted = verify_decomposition(decomposition, output=z, domain=sp.And(sp.Ne(x, 0), y >= 0))
    assert restricted.equivalent and restricted.certified and restricted.function_verified


@pytest.mark.parametrize(
    "kwargs",
    [
        {"output": x},
        {"output": t},
        {"outer_selector": x > 0},
        {"inner_selector": z > 0},
        {"target_selector": t > 0},
        {"domain": z > 0},
    ],
)
def test_selector_scope(kwargs):
    decomposition = RationalDecomposition(
        (x * x + y * y) ** 2, (x, y), t, sp.QQ, t * t, x * x + y * y
    )
    options = {"output": z}
    options.update(kwargs)
    with pytest.raises(ValueError):
        verify_decomposition(decomposition, **options)


def test_invalid_certificate():
    decomposition = RationalDecomposition(
        (x * x + y * y) ** 2, (x, y), t, sp.QQ, t * t, x * x + y * y
    )
    with pytest.raises(ValueError, match="failed verification"):
        verify_decomposition(replace(decomposition, outer=t * t + 1), output=z)
    constants = sp.QQ.algebraic_field(sp.I)
    complex_certificate = RationalDecomposition(
        (x * x + y * y) ** 2 + sp.I, (x, y), t, constants, t * t + sp.I, x * x + y * y
    )
    with pytest.raises(ValueError, match="real embedding"):
        verify_decomposition(complex_certificate, output=z)
