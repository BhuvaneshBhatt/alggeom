"""Exact recursive pointwise covers over nonprime parameter quotients."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.localized_linear_cover import (
    LinearCoverLimitError,
    comprehensive_localized_linear_cover,
)
from alggeom.parameter_localization import ParameterLocalization


def test_exceptional_pivot():
    a, b, x = sp.symbols("a b x")
    ring = PolynomialRing((a, b))
    cover = comprehensive_localized_linear_cover(
        (a * x + b, b * x + a), x, ParameterLocalization(ring.ideal())
    )
    assert len(cover.pieces) == 4
    assert cover.verify()
    assert cover.specialize((1, 1)) == (x + 1,)
    assert cover.specialize((2, 0)) == (1,)
    assert cover.specialize((0, 0)) == ()
    assert not replace(cover, pieces=cover.pieces[:-1]).verify()
    with pytest.raises(LinearCoverLimitError):
        comprehensive_localized_linear_cover(
            (a * x + b, b * x + a), x, ParameterLocalization(ring.ideal()), max_pieces=2
        )


@pytest.mark.parametrize("relations", [("a*b",), ("a*b", "b**2")])
def test_nonreduced_parameters(relations):
    a, b, x = sp.symbols("a b x")
    ring = PolynomialRing((a, b))
    ideal = ring.ideal(*(sp.sympify(value, locals={"a": a, "b": b}) for value in relations))
    cover = comprehensive_localized_linear_cover(
        (a * x + b, b * x + a), x, ParameterLocalization(ideal)
    )
    assert cover.verify()
    assert cover.specialize((0, 0)) == ()
    assert cover.specialize((2, 0)) == (1,)
    assert not replace(cover, pieces=cover.pieces + cover.pieces[:1]).verify()


def test_localized_specialization():
    a, b, x = sp.symbols("a b x")
    ring = PolynomialRing((a, b))
    cover = comprehensive_localized_linear_cover(
        (a * x + b, (a**2 - 2) * x), x, ParameterLocalization(ring.ideal(), (a + b,))
    )
    assert cover.verify()
    assert sp.simplify(cover.specialize((sp.sqrt(2), 1))[0] - (x + sp.sqrt(2) / 2)) == 0
    with pytest.raises(ValueError, match="not exactly one"):
        cover.specialize((0, 0))


def test_zero_generators_and_non_linear_failure():
    a, x = sp.symbols("a x")
    base = ParameterLocalization(PolynomialRing((a,)).ideal(a * a))
    cover = comprehensive_localized_linear_cover((a * x,), x, base)
    assert cover.verify()
    assert cover.specialize((0,)) == ()
    with pytest.raises(NotImplementedError, match="degree"):
        comprehensive_localized_linear_cover((a * x * x,), x, base)


def test_ordered_exceptional_ideals():
    a, b, c, x = sp.symbols("a b c x")
    ring = PolynomialRing((a, b, c))
    cover = comprehensive_localized_linear_cover(
        (a * x + 1, b * x + 1, c * x + 1), x, ParameterLocalization(ring.ideal(a * b))
    )
    assert cover.verify()
    assert cover.specialize((2, 0, 2)) == (1,)
    assert cover.specialize((2, 0, 0)) == (1,)
    assert cover.specialize((0, 0, 0)) == (1,)
