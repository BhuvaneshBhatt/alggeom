"""Independent mathematical answers for the application workflows."""

from pathlib import Path

import pytest
import sympy as sp

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def test_equation_roots(run_example):
    result = run_example(EXAMPLES / "solve_equation.py")
    expected_real = {sp.S.Zero, -sp.S.One, (-1 + sp.sqrt(5)) / 2, (-1 - sp.sqrt(5)) / 2}
    expected = expected_real | {(-1 + sp.sqrt(3) * sp.I) / 2, (-1 - sp.sqrt(3) * sp.I) / 2}
    assert result["solutions"] == expected
    assert result["real_solutions"] == expected_real
    assert sorted(sp.degree(f.expr, f.variables[0]) for f in result["components"]) == [2, 3]


@pytest.mark.real
def test_response_contour(run_example):
    result = run_example(EXAMPLES / "response_contours.py")
    x, y = result["x"], result["y"]
    radius = x * x + y * y
    numerator, denominator = sp.cancel(2 * result["response"] - 1).as_numer_denom()
    assert sp.expand(numerator - (radius - 1) * (radius + 1)) == 0
    assert sp.expand(denominator - (1 + radius**2)) == 0
    assert len(result["scalar_candidates"]) == 2
    assert len(result["admissible_levels"]) == 1
    parameter = sp.Symbol("parameter", real=True)
    circle = {x: (1 - parameter**2) / (1 + parameter**2), y: 2 * parameter / (1 + parameter**2)}
    assert sp.cancel(result["response"].subs(circle) - sp.Rational(1, 2)) == 0
    assert result["response"].subs({x: 0, y: 0}) == 0


def test_curve_inverse(run_example):
    result = run_example(EXAMPLES / "curve_reparametrization.py")
    t, s, u, v = (result[name] for name in ("t", "s", "u", "v"))
    assert all(
        coordinate.subs(t, -t) == coordinate for coordinate in result["original_coordinates"]
    )
    assert not result["field"].contains(t)
    # Check the inverse on the implicit curve, not just on its parametrized image.
    curve = sp.groebner([v * v - u**3], v, u, domain=sp.QQ)
    for coordinate, expected in zip(result["new_coordinates"], (u, v), strict=True):
        numerator = sp.cancel(
            coordinate.subs(s, result["inverse_parameter"]) - expected
        ).as_numer_denom()[0]
        assert curve.reduce(numerator)[1] == 0


@pytest.mark.real
@pytest.mark.parametrize("a,b", [(-2, 3), (0, -1), (2, 3)])
@pytest.mark.parametrize("offset", [0, 1])
def test_quadratic_conditions(run_example, a, b, offset):
    result = run_example(EXAMPLES / "quadratic_towers.py")
    x = result["x"]
    c = sp.Rational(a * b, 2) - sp.Rational(a**3, 8) + offset
    assignment = {result["a"]: a, result["b"]: b, result["c"]: c}
    specialized = sp.expand(result["polynomial"].subs(assignment))
    reference = sp.decompose(specialized, x)
    assert (result["pair"].condition.subs(assignment) is sp.true) == (len(reference) > 1)
    assert len(reference) == (2 if offset == 0 else 1)
