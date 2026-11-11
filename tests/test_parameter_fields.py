"""Exact quotient-field arithmetic and recursive cover certification."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.parameter_cover import comprehensive_prime_cover
from alggeom.parameter_fields import (
    ParameterFieldError,
    PrimeParameterField,
    groebner_over_prime_parameters,
)


def test_implicit_fraction_field():
    a, b = sp.symbols("a b")
    field = PrimeParameterField(PolynomialRing((a, b)).ideal(a**2 + b**2))
    assert field.reduce(a**2 + b**2) == 0
    assert field.element(a) ** 2 == field.element(-(b**2))
    assert field.element(a) ** 2 / field.element(b) == field.element(-b)
    assert field.element(a).inverse() * field.element(a) == field.element(1)
    with pytest.raises(ZeroDivisionError):
        field.element(1, a**2 + b**2)


def test_nonprime_rejected():
    a, b = sp.symbols("a b")
    for equations in ((a**2,), (a * b,)):
        with pytest.raises(ParameterFieldError):
            PrimeParameterField(PolynomialRing((a, b)).ideal(*equations))


def test_quotient_basis_certificate():
    a, b, x, y = sp.symbols("a b x y")
    prime = PolynomialRing((a, b)).ideal(a**2 + b**2)
    result = groebner_over_prime_parameters((a * x + b * y, y**2 - 1), (x, y), prime)
    assert result.verify()
    assert a in result.denominator_guards
    assert not replace(result, denominator_guards=()).verify()
    corrupted = replace(
        result, membership_witnesses=(((sp.S.Zero,) * 2), *result.membership_witnesses[1:])
    )
    assert not corrupted.verify()


def test_recursive_implicit_cover():
    a, b, x, y = sp.symbols("a b x y")
    equations = ((a**2 + b**2 + 1) * x, a * y)
    cover = comprehensive_prime_cover(equations, (x, y), (a, b))
    assert len(cover.pieces) == 4
    assert cover.verify()
    for values in ((1, 1), (0, 1), (0, sp.I), (sp.I, 0), (2, 3)):
        substitutions = dict(zip((a, b), values, strict=True))
        expected = sp.groebner(
            tuple(expr.subs(substitutions) for expr in equations),
            x,
            y,
            order="lex",
            extension=True,
        )
        actual = sp.groebner(cover.specialize(values), x, y, order="lex", extension=True)
        assert tuple(poly.as_expr() for poly in actual.polys) == tuple(
            poly.as_expr() for poly in expected.polys
        )


def test_certificate_rejects_wrong_stratum():
    a, b, x = sp.symbols("a b x")
    cover = comprehensive_prime_cover(((a**2 + b**2) * x,), (x,), (a, b))
    assert cover.verify()
    branch = next(
        piece for piece in cover.pieces if piece.certificate.parameter_ideal.is_zero is False
    )
    modified = replace(branch, condition=sp.true)
    assert not modified.verify()


def test_recursive_budget():
    a, b, x, y = sp.symbols("a b x y")
    from alggeom.parameter_cover import ParameterCoverLimitError

    with pytest.raises(ParameterCoverLimitError):
        comprehensive_prime_cover((a * x, b * y), (x, y), (a, b), max_nodes=1)


def test_nonprincipal_parameter_prime():
    a, b, c, x, y = sp.symbols("a b c x y")
    prime = PolynomialRing((a, b, c)).ideal(c, a**2 + b**2)
    field = PrimeParameterField(prime)
    assert field.element(a**2) == field.element(-(b**2))
    assert field.reduce(c) == 0
    proof = groebner_over_prime_parameters((a * x + b * y, y**2 - 1), (x, y), prime)
    assert proof.verify()
    assert a in proof.denominator_guards


def test_nested_nonprincipal_branch():
    a, b, c, x, y = sp.symbols("a b c x y")
    equations = ((a**2 + b**2) * x, c * y)
    cover = comprehensive_prime_cover(equations, (x, y), (a, b, c))
    assert cover.verify()
    assert len(cover.pieces) == 4
    for values, target in (
        ((1, 1, 1), (x, y)),
        ((1, 1, 0), (x,)),
        ((sp.I, 1, 1), (y,)),
        ((sp.I, 1, 0), ()),
    ):
        selected = cover.specialize(values)
        assert tuple(sp.Poly(item, x, y).monic().as_expr() for item in selected) == target
