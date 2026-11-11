"""Exact constructible projection, including nonlinear and excluded-locus cases."""

import pytest
import sympy as sp

from alggeom.chevalley import ProjectionLimitError, project_constructible
from alggeom.constructible import ConstructibleSet, LocallyClosedSet
from alggeom.ideals import PolynomialRing


def project(equations, variables, eliminate, *, nonzero=()):
    ring = PolynomialRing(tuple(variables))
    source = ConstructibleSet((LocallyClosedSet(ring.ideal(*equations), nonzero),))
    return project_constructible(source, eliminate)


def test_nonclosed_polynomial_image():
    s, t, x, y = sp.symbols("s t x y")
    image = project((x - s, y - s * t), (s, t, x, y), (s, t))
    assert image.contains((0, 0))
    assert image.contains((2, 5))
    assert not image.contains((0, 1))
    assert image.closure_ideal.is_zero


def test_nonlinear_image():
    s, t, x, y = sp.symbols("s t x y")
    image = project((x - s**2, y - s * t), (s, t, x, y), (s, t))
    assert image.contains((0, 0))
    assert image.contains((1, 3))
    assert image.contains((4, -2))
    assert not image.contains((0, 3))


def test_repeated_root_and_removed_origin():
    a, x = sp.symbols("a x")
    image = project((x**2 - a,), (x, a), (x,), nonzero=(x,))
    assert image.contains((2,))
    assert image.contains((-2,))
    assert not image.contains((0,))


def test_parameter_degree_drop():
    a, x = sp.symbols("a x")
    image = project((a * x - 1,), (x, a), (x,))
    assert not image.contains((0,))
    assert image.contains((1,))


def test_simultaneous_nonlinear_equations():
    a, b, x = sp.symbols("a b x")
    image = project((x**2 - a, x**2 - b), (x, a, b), (x,))
    assert image.contains((2, 2))
    assert image.contains((0, 0))
    assert not image.contains((2, 3))


def test_multiple_source_variables():
    a, x, y = sp.symbols("a x y")
    image = project((x**2 + y**2 - a,), (x, y, a), (x, y))
    assert image.contains((0,))
    assert image.contains((3,))
    assert image.contains((-1,))


def test_budget_is_explicit():
    a, x = sp.symbols("a x")
    ring = PolynomialRing((x, a))
    source = ConstructibleSet((LocallyClosedSet(ring.ideal(x**2 - a)),))
    with pytest.raises(ProjectionLimitError):
        project_constructible(source, (x,), max_branches=1)


def test_parameterized_gcd_specializations():
    """Cross-check constructible strata against independent specialized gcds."""
    from functools import reduce

    a, b, x = sp.symbols("a b x")
    ring = PolynomialRing((x, a, b))
    examples = (
        ((a * x * x + b * x + 1, x * x - a), x - b),
        ((a * x * x, b * x * x + 1), x),
        ((x * x - a, x * x - b), x * x),
        ((a * x * x + b * x, x * x - b), x - a),
        ((a * x * x + b, x + a), x),
    )
    for equations, guard in examples:
        pieces = project_constructible(
            ConstructibleSet((LocallyClosedSet(ring.ideal(*equations), (guard,)),)),
            (x,),
        ).pieces
        for a_value in (-1, 0, 1):
            for b_value in (-1, 0, 1):
                substitutions = {a: a_value, b: b_value}
                specialized = [sp.Poly(f.subs(substitutions), x, domain=sp.QQ) for f in equations]
                if all(f.is_zero for f in specialized):
                    expected = not sp.Poly(guard.subs(substitutions), x).is_zero
                else:
                    gcd = reduce(sp.gcd, specialized)
                    excluded = sp.Poly(guard.subs(substitutions), x, domain=sp.QQ)
                    expected = (
                        gcd.degree() > 0
                        and sp.gcd(gcd, excluded ** gcd.degree()).degree() < gcd.degree()
                    )
                actual = any(
                    all(f.subs(substitutions) == 0 for f in piece.ideal.generators)
                    and all(g.subs(substitutions) != 0 for g in piece.nonzero)
                    for piece in pieces
                )
                assert actual == expected, (equations, guard, substitutions)
