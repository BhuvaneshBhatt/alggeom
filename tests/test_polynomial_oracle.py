"""Compare complete polynomial chains with SymPy's independent polynomial algorithm."""

import pytest
import sympy as sp

from alggeom.decomposition import decompose

x = sp.Symbol("x")
quadratic = x**2 + 2 * x + 3
POLYNOMIALS = [
    x**6,
    x**4 + x,
    x**6 + x,
    sp.chebyshevt(6, x),
    quadratic**3 - 3 * quadratic + 1,
    ((x**2 + x) ** 2 + 1) ** 2,
    (x**3 + x) ** 2 + 2 * (x**3 + x),
    sp.Rational(1, 2) * (x**2 + x / 3) ** 3 + x**2 + x / 3,
]


def _map_degree(component):
    numerator, denominator = sp.fraction(sp.cancel(component.expr))
    variable = component.variables[0]
    return max(sp.degree(numerator, variable), sp.degree(denominator, variable))


@pytest.mark.parametrize(
    "polynomial",
    POLYNOMIALS,
    ids=[
        "power",
        "quartic",
        "sextic",
        "chebyshev",
        "shifted",
        "tower",
        "cubic_inner",
        "rational_coefficients",
    ],
)
def test_polynomial_chain(polynomial):
    polynomial = sp.expand(polynomial)
    reference = sp.decompose(polynomial, x)
    actual = decompose(polynomial, x)
    reference_value = x
    for component in reversed(reference):
        reference_value = component.subs(x, reference_value)
    actual_value = x
    for component in reversed(actual):
        actual_value = component(actual_value)
    assert sp.expand(reference_value - polynomial) == 0
    assert sp.cancel(actual_value - polynomial) == 0
    # Complete polynomial chains in characteristic zero have the same degree multiset.
    assert sorted(_map_degree(component) for component in actual) == sorted(
        sp.degree(component, x) for component in reference
    )
