import pytest
import sympy as sp
from algroots.domains import field_element
from algroots.finite_decomposition import nilradical_basis, primary_decomposition
from algroots.finite_fields import factor_univariate
from algroots.quotient import QuotientAlgebra

from alggeom.fields import field_order, finite_field
from alggeom.ideals import PolynomialRing
from alggeom.resolutions import free_resolution
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, a = sp.symbols("x y a")


def test_extension_field():
    field = finite_field(3, 2)
    assert field_order(field) == 9
    generator = field.generator
    assert generator**9 == generator
    assert generator / generator == field.one
    assert field.exquo(field.one, generator) * generator == field.one
    polynomial = sp.Poly(x * x + 1, x, domain=field)
    unit, factors = factor_univariate(polynomial)
    product = sp.Poly(unit, x, domain=field)
    for factor, exponent in factors:
        product *= factor**exponent
    assert product == polynomial
    assert [f.degree() for f, e in factors] == [1, 1]
    with pytest.raises(ValueError, match="irreducible"):
        finite_field(3, modulus=a * a - 1)
    with pytest.raises(ValueError, match="prime"):
        finite_field(9)
    with pytest.raises(ValueError, match="distinct"):
        PolynomialRing((field.symbol,), field)


def test_frobenius_nilpotents():
    algebra = QuotientAlgebra.from_polynomials((x * x, y * y), (x, y), domain=sp.GF(2))
    assert algebra.trace_pairing == sp.zeros(4)
    assert algebra.geometric_solution_count == 1
    assert len(nilradical_basis(algebra)) == 3
    assert not algebra.is_field
    split = QuotientAlgebra.from_polynomials((x * x - x, y * y - y), (x, y), domain=sp.GF(2))
    factors = primary_decomposition(split)
    assert len(factors) == 4
    assert all(c.dimension == r.dimension == 1 and r.is_field for c, r in factors)
    extension = QuotientAlgebra.from_polynomials((x * x + x + 1,), (x,), domain=sp.GF(2))
    assert extension.is_field and extension.geometric_solution_count == 2


@pytest.mark.parametrize("field", [sp.GF(2), finite_field(3, 2)])
def test_finite_geometry(field):
    q = field_order(field)
    ring = PolynomialRing((x, y), field)
    assert len(AffineScheme(ring).rational_points()) == q * q
    assert len(ProjectiveScheme(ring).rational_points()) == q + 1
    ideal = ring.ideal((x * y) ** 2)
    assert ideal.radical() == ring.ideal(x * y)
    fat = ring.ideal(x**4, y * y)
    assert fat.radical() == ring.ideal(x, y)
    assert len(fat.primary_decomposition()) == 1
    resolution = free_resolution(ring.ideal(x, y))
    assert resolution.shifts == ((0,), (1, 1), (2,)) and resolution.verify()
    with pytest.raises(NotImplementedError, match="positive-dimensional"):
        ideal.primary_decomposition()
    with pytest.raises(ValueError, match="characteristic"):
        PolynomialRing((x,)).quotient().base_change(field)


def test_height_search():
    curve = AffineScheme(PolynomialRing((x, y)).ideal(y - x * x))
    assert {p.coordinates for p in curve.point_search(2)} == {(-1, 1), (0, 0), (1, 1)}
    line = ProjectiveScheme(PolynomialRing((x, y)))
    assert {p.coordinates for p in line.point_search(1)} == {(1, -1), (1, 0), (1, 1), (0, 1)}
    assert len(line.point_search(2)) == 8
    with pytest.raises(ValueError, match="height"):
        curve.point_search()
    with pytest.raises(RuntimeError, match="max_candidates"):
        curve.point_search(100, max_candidates=10)
    with pytest.raises(RuntimeError, match="max_candidates"):
        AffineScheme(PolynomialRing((x, y), sp.GF(101))).rational_points(max_candidates=100)


def test_finite_coefficients():
    field = finite_field(5, 2)
    a = field.symbol
    coefficient = (a + 1) / 2
    value = field_element(coefficient, field)
    assert value * field.convert(2) == field.convert(a + 1)
    ring = PolynomialRing((x,), field)
    algebra = ring.quotient(x * x - coefficient).finite_algebra()
    parameter = sp.Symbol("T")
    assert algebra.minimal_polynomial(x, parameter) == sp.Poly(
        parameter**2 - field.to_sympy(value), parameter, domain=field
    )
    assert algebra.trace(x) == 0
    assert algebra.norm(x) == field.to_sympy(-value)
    assert algebra.multiplication_matrix(x**4) == sp.eye(2) * field.to_sympy(value**2)
    generic = PolynomialRing((x,), field.frac_field(y)).quotient(
        x * x - y - coefficient, groebner_backend="sympy"
    )
    assert generic(x * x) == generic(field.to_sympy(value) + y)


def test_field_embeddings():
    field = finite_field(2, 2)
    source = PolynomialRing((x,), sp.GF(2)).quotient(x * x + x + 1)
    extended = source.base_change(field)
    assert len(primary_decomposition(extended.finite_algebra())) == 2
