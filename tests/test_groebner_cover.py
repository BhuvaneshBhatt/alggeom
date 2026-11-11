"""Gröbner bases must remain valid on every claimed specialization stratum."""

import sympy as sp

from alggeom.groebner_cover import (
    comprehensive_groebner_cover,
)


def test_linear_specialization():
    a, x, y = sp.symbols("a x y")
    cover = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,))
    assert cover.verify_at((0,))
    assert cover.verify_at((1,))
    assert cover.verify_at((-3,))
    assert cover.specialize((0,)) == (x + y,)


def test_two_parameter_strata():
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover((a * x, b * y), (x, y), (a, b))
    assert len(cover.pieces) == 4
    for values in ((0, 0), (0, 1), (2, 0), (2, 3), (-2, 0)):
        assert cover.verify_at(values)


def test_coupled_nonlinear_family():
    a, x, y = sp.symbols("a x y")
    cover = comprehensive_groebner_cover((x * y - a, y - x), (x, y), (a,))
    for value in (0, 1, 2, -1):
        assert cover.verify_at((value,))


def test_affine_parameter_diagonal():
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover(((a + b) * x, x * y - 1), (x, y), (a, b))
    for value in ((0, 0), (1, -1), (1, 2), (2, -2)):
        assert cover.verify_at(value)


def test_nonlinear_exceptional_hypersurface():
    a, b, x = sp.symbols("a b x")
    cover = comprehensive_groebner_cover(((a * a + b * b) * x,), (x,), (a, b))
    assert len(cover.pieces) >= 2
    for values in ((0, 0), (1, 0), (0, 2), (2, 3), (sp.I, 1)):
        assert cover.verify_at(values)


def test_nested_nonlinear_branch():
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover(((a * a + b * b) * x, a * y), (x, y), (a, b))
    for values in ((0, 0), (0, 2), (2, 0), (sp.I, 1), (-sp.I, 1)):
        assert cover.verify_at(values)
