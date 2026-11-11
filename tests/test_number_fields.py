"""Exact embeddings, field-dependent answers and reconstruction over number fields."""

import pytest
import sympy as sp
from algroots.function_fields import RationalCoordinateMap
from algroots.quotient import QuotientAlgebra
from algroots.subfields import field_normal_form

from alggeom.decomposition.algebraic import reconstruct_algebraic
from alggeom.decomposition.factors import discover_factors
from alggeom.decomposition.rational import RationalDecomposition, search_rational

x, z, a, t = sp.symbols("x z a t")
REAL_GENERATORS = [sp.sqrt(3), sp.real_root(2, 3), sp.sqrt(2) + sp.sqrt(3)]


@pytest.mark.parametrize("c", [2, 3, 5, -3])
def test_quadratic_splitting(c):
    alpha = sp.sqrt(c)
    constants = sp.QQ.algebraic_field(alpha)
    function = x**2 + c / x**2
    rational = search_rational(function, x, parameter=t)
    extended = search_rational(function, x, parameter=t, constants=constants)
    assert rational.complete and rational.verify() and len(rational.decompositions) == 1
    assert extended.complete and extended.verify() and len(extended.decompositions) == 3
    expected = [x**2, x + alpha / x, x - alpha / x]
    ambient = constants.frac_field(x)
    for generator in expected:
        matches = []
        for result in extended.decompositions:
            forward = RationalCoordinateMap(ambient, generator, t).express(result.inner)
            inverse = RationalCoordinateMap(ambient, result.inner, t).express(generator)
            if forward is not None and inverse is not None:
                matches.append(result)
        assert len(matches) == 1
    assert rational.decompositions[0].inner == x**2


@pytest.mark.parametrize("alpha", REAL_GENERATORS, ids=["quadratic", "cubic", "composite"])
def test_embedded_coordinates(alpha):
    constants = sp.QQ.algebraic_field(alpha)
    inner = x**2 + alpha * x
    certificate = RationalDecomposition(inner**2 + alpha, (x,), t, constants, t**2 + alpha, inner)
    assert certificate.verify()
    coordinates = RationalCoordinateMap(constants.frac_field(x), inner, t)
    assert coordinates.express(certificate.function) == t**2 + alpha
    assert coordinates.express(x) is None
    result = search_rational(x**4 + alpha, x, parameter=t, constants=constants)
    assert result.complete and result.verify()
    assert [(d.outer, d.inner) for d in result.decompositions] == [(t**2 + alpha, x**2)]


@pytest.mark.parametrize("alpha", REAL_GENERATORS, ids=["quadratic", "cubic", "composite"])
def test_implicit_fields(alpha):
    constants = sp.QQ.algebraic_field(alpha)
    ambient = QuotientAlgebra.from_polynomials([a - z], (a,), domain=constants.frac_field(z))
    relation = x**2 - z**4 - alpha
    result = reconstruct_algebraic(
        relation, (x, z), candidate=relation, ambient=ambient, parameter=t, constants=constants
    )
    assert result.verify() and result.status == "decomposable"
    decomposition = result.decomposition
    assert decomposition.constants == constants and decomposition.exponent == 1
    resultant = sp.Poly(
        sp.resultant(decomposition.outer.as_expr(), decomposition.inner.as_expr(), t),
        x,
        z,
        domain=constants,
    )
    assert resultant == decomposition.unit * sp.Poly(relation, x, z, domain=constants)


@pytest.mark.parametrize("alpha", REAL_GENERATORS, ids=["quadratic", "cubic", "composite"])
def test_factor_extensions(alpha):
    constants = sp.QQ.algebraic_field(alpha)
    relation = x**2 + alpha * z * x + z
    discovery = discover_factors(relation, (x, z), constants=constants)
    assert discovery.complete and discovery.verify()
    assert [orbit.extension_degree for orbit in discovery.orbits] == [2, 1]
    for orbit in discovery.orbits:
        assert orbit.ambient.domain == constants.frac_field(z)
        difference = sp.Poly(orbit.factor * orbit.cofactor - relation, x)
        assert all(field_normal_form(orbit.ambient, c) == 0 for c in difference.all_coeffs())


def test_imaginary_quotient():
    constants = sp.QQ.algebraic_field(sp.sqrt(-3))
    ambient = QuotientAlgebra.from_polynomials([a - z], (a,), domain=constants.frac_field(z))
    relation = x**2 - z**4 - sp.sqrt(-3)
    result = reconstruct_algebraic(
        relation, (x, z), candidate=relation, ambient=ambient, parameter=t, constants=constants
    )
    assert result.verify() and result.decomposition.constants == constants


@pytest.mark.parametrize(
    "alpha", [sp.sqrt(-3), sp.sqrt(2), sp.real_root(2, 3), sp.sqrt(2) + sp.sqrt(3)]
)
def test_quotient_arithmetic(alpha):
    constants = sp.QQ.algebraic_field(alpha)
    base = constants.frac_field(z)
    ambient = QuotientAlgebra.from_polynomials([a**2 - z - alpha], (a,), domain=base)
    assert ambient.is_field and ambient.dimension == 2
    minimum = ambient.minimal_polynomial(a, t)
    assert (
        minimum.domain == base
        and minimum.gens == (t,)
        and minimum.all_coeffs() == [1, 0, -z - alpha]
    )
    assert ambient.trace(a) == 0 and ambient.normal_form(ambient.norm(a) + z + alpha) == 0
    assert ambient.normal_form(a * ambient.inverse(a)) == 1
    reduced = ambient.normal_form((a**2 + alpha * a) / (z + 1))
    assert sp.cancel(reduced - (z + alpha + alpha * a) / (z + 1)) == 0


def test_rootof_bound_coordinate():
    from alggeom.ideals import PolynomialRing

    x = sp.Symbol("x")
    alpha = sp.CRootOf(x**3 - x - 1, 0)
    field = sp.QQ.algebraic_field(alpha)
    ring = PolynomialRing((x,), field)
    ideal = ring.ideal(x - alpha)
    assert ideal.contains(x - alpha)
    assert ideal.contains(x**2 - 2 * alpha * x + alpha**2)
    assert not ideal.contains(x)


def test_rootof_rejects_inexact_coefficient():
    from alggeom.ideals import PolynomialRing

    x = sp.Symbol("x")
    alpha = sp.CRootOf(x**3 - x - 1, 0)
    ring = PolynomialRing((x,), sp.QQ.algebraic_field(alpha))
    with pytest.raises(ValueError, match="exact"):
        ring.ideal(x - alpha + sp.Float("0.125"))
