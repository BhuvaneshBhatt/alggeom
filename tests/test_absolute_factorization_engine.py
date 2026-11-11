"""Independent reconstruction and geometric-component regression cases."""

import pytest
import sympy as sp

from alggeom.absolute_factorization import absolute_factor_bivariate
from alggeom.geometric import absolute_primary_decomposition, is_geometrically_integral
from alggeom.ideals import PolynomialRing


@pytest.mark.parametrize(
    ("expression", "factor_count"),
    [
        (lambda x, y: x**2 + y**2, 2),
        (lambda x, y: x**3 + y**3, 3),
        (lambda x, y: x**4 + y**4, 4),
        (lambda x, y: x**2 + y**2 + 1, 1),
        (lambda x, y: (x**2 + y**2) ** 2, 2),
        (lambda x, y: (x**2 + y**2) * (x + y + 1), 3),
        (lambda x, y: x * y * (x + y), 3),
    ],
)
def test_absolute_factorization(expression, factor_count):
    x, y = sp.symbols("x y")
    polynomial = expression(x, y)
    factors = absolute_factor_bivariate(polynomial, (x, y))
    assert len(factors.factors) == factor_count
    recovered = sp.Poly(
        sp.prod(base**power for base, power in factors.factors),
        x,
        y,
        domain=factors.field,
    )
    expected = sp.Poly(polynomial, x, y, domain=factors.field)
    assert expected.exquo(recovered).total_degree() == 0
    for index, (factor, _) in enumerate(factors.factors):
        for other, _ in factors.factors[index + 1 :]:
            assert (
                sp.gcd(
                    sp.Poly(factor, x, y, domain=factors.field),
                    sp.Poly(other, x, y, domain=factors.field),
                ).total_degree()
                == 0
            )


def test_primary_reconstruction_with_multiplicity():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    ideal = ring.ideal((x**2 + y**2) ** 2 * (x + y + 1))
    result = absolute_primary_decomposition(ideal)
    assert len(result.components) == 3
    recovered = result.components[0].ring.ideal(1)
    for component in result.components:
        recovered = recovered.intersection(component)
    assert recovered == ideal.change_ring(result.components[0].ring)
    assert not is_geometrically_integral(ring.ideal(x**4 + y**4))
    assert is_geometrically_integral(ring.ideal(y**2 - x**3))


def test_three_variable_geometric_primary_factors():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal((x**2 + y**2) ** 2 * (z + 1))
    result = absolute_primary_decomposition(ideal)
    assert len(result.components) == 3
    recovered = result.components[0].ring.ideal(1)
    for component in result.components:
        recovered = recovered.intersection(component)
    assert recovered == ideal.change_ring(result.components[0].ring)
    assert is_geometrically_integral(ring.ideal(x**2 + y**2 + z**2))
