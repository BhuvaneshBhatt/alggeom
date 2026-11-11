"""Fixed-degree polynomial conditions against independent coefficient elimination."""

import pytest
import sympy as sp
from semialg import equivalent, quantifier_eliminate

from alggeom.decomposition.parameters import polynomial_pair, polynomial_system

x, y, z, a, b, c, t = sp.symbols("x y z a b c t", real=True)


def test_quartic_condition():
    result = polynomial_pair(
        x**4 + a * x**3 + b * x * x + c * x,
        x,
        degrees=(2, 2),
        parameters=(a, b, c),
        component_variable=t,
    )
    expected = sp.Eq(a**3 - 4 * a * b + 8 * c, 0)
    assert equivalent(result.condition, expected, (a, b, c))
    system = result.system
    projected = quantifier_eliminate(
        system.formula,
        [("exists", v) for v in system.unknowns],
        variables=(*system.parameters, *system.unknowns),
    )
    assert equivalent(result.condition, projected, (a, b, c))
    assert result.inner == x * x + a * x / 2
    assert sp.expand(result.outer - (t * t + (b - a * a / 4) * t)) == 0
    assert (
        sp.expand((result.composition - system.polynomial.as_expr()).subs(c, a * b / 2 - a**3 / 8))
        == 0
    )
    assert equivalent(result.condition.subs(a, 0), sp.Eq(c, 0), (b, c))
    assignment = dict(result.coefficient_values)
    assert set(assignment) == set(system.unknowns)
    assert sp.simplify(system.formula.subs(assignment).subs(c, a * b / 2 - a**3 / 8)) is sp.true


@pytest.mark.parametrize("degrees", [(2, 3), (3, 2), (3, 3), (2, 4)])
def test_normalized_components(degrees):
    m, n = degrees
    inner = 2 * x**n + x ** (n - 1) - 3 * x + 5
    outer = 3 * t**m - 2 * t + 7
    f = sp.expand(outer.subs(t, inner))
    result = polynomial_pair(f, x, degrees=degrees, component_variable=t)
    assert result.condition is sp.true
    assert sp.expand(result.composition - f) == 0
    assert sp.Poly(result.inner, x).LC() == 1
    assert result.inner.subs(x, 0) == 0
    assert sp.degree(result.outer, t) == m


def test_negative_condition():
    result = polynomial_pair(x**4 + x, x, degrees=(2, 2))
    assert result.condition is sp.false


def test_degree_and_pole_guards():
    result = polynomial_pair(a * x**4 + b * x * x + 1, x, degrees=(2, 2), parameters=(a, b))
    assert equivalent(result.condition, sp.Ne(a, 0), (a, b))
    rational = polynomial_pair((x * x + 1) ** 2 / a, x, degrees=(2, 2), parameters=(a,))
    assert rational.condition.subs(a, 0) is sp.false
    assert rational.condition.subs(a, 2) is sp.true
    assert sp.cancel(rational.composition - (x * x + 1) ** 2 / a) == 0


def test_algebraic_constants():
    f = sp.expand((x * x + sp.sqrt(2) * x) ** 2 + 3 * (x * x + sp.sqrt(2) * x))
    result = polynomial_pair(f, x, degrees=(2, 2))
    assert result.condition is sp.true
    assert sp.expand(result.composition - f) == 0


@pytest.mark.parametrize(
    "f,degrees,params",
    [
        (x**3, (2, 2), ()),
        (x**4, (1, 4), ()),
        (x**4, (True, 4), ()),
        (x**4 + 0.5, (2, 2), ()),
        (x**4 + a, (2, 2), ()),
        (1 / x, (2, 2), ()),
        (x**4 + sp.I, (2, 2), ()),
        (x**4 + sp.pi, (2, 2), ()),
        (x**4, (2, 2), (x,)),
        (x**4, (2, 2), (a, a)),
    ],
)
def test_invalid_requests(f, degrees, params):
    with pytest.raises(ValueError):
        polynomial_system(f, x, degrees=degrees, parameters=params)
