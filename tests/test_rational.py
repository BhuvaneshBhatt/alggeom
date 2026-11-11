"""Rational two-map classes, field equivalence, candidate coverage and certificates."""

from dataclasses import replace
from itertools import combinations

import pytest
import sympy as sp
from algroots.function_fields import luroth_generator

from alggeom.decomposition import DecompositionError
from alggeom.decomposition.rational import RationalDecomposition, RationalSearch, search_rational

x, y, z, t = sp.symbols("x y z t")


def assert_compositions(result):
    assert result.verify()
    for decomposition in result.decompositions:
        assert decomposition.verify()
        assert (
            sp.cancel(decomposition.outer.subs(t, decomposition.inner) - result.original_expression)
            == 0
        )


@pytest.mark.parametrize(
    "function,degrees",
    [
        ((x * x + y * y) ** 2, {(2, 2)}),
        ((x * y + z) ** 3, {(3, 2)}),
        ((x + y) ** 6, {(3, 2), (2, 3)}),
        (
            ((x * x + y * y) ** 2 + (x * y + 1) ** 2)
            / ((x * x + y * y) ** 2 + (x * x + y * y) * (x * y + 1) + 2 * (x * y + 1) ** 2),
            {(2, 2)},
        ),
    ],
)
def test_multivariate_pairs(function, degrees):
    result = search_rational(function, (x, y, z), parameter=t)
    assert isinstance(result, RationalSearch)
    assert result.complete and result.status == "decomposable"
    assert {d.degrees for d in result.decompositions} == degrees
    assert_compositions(result)
    assert not replace(result, decompositions=()).verify()


def test_nonlinear_diagonal():
    result = search_rational((x * x + y * y) ** 2, (x, y), parameter=t)
    assert result.candidate_count == 1 and len(result.factors) == 2
    assert all(p.total_degree() == 2 for p in result.factors)
    assert result.decompositions[0].inner == x * x + y * y
    assert_compositions(result)


def test_published_example():
    a, b = (x * y + x * z * z - y**5, 9 * x + 25 * y * z)
    function = (a * a - b * b) / (a * a - 2 * b * b)
    result = search_rational(function, (x, y, z), parameter=t)
    assert result.complete and result.status == "decomposable"
    assert result.candidate_count == 1
    assert [d.degrees for d in result.decompositions] == [(2, 5)]
    assert_compositions(result)
    field = luroth_generator([result.decompositions[0].inner], variables=(x, y, z))
    assert field.contains(a / b)


def test_distinct_fields():
    result = search_rational((x / y + y / x) ** 2, (x, y), parameter=t)
    assert result.complete and len(result.decompositions) == 3
    assert {d.degrees for d in result.decompositions} == {(2, 2)}
    assert_compositions(result)
    fields = [luroth_generator([d.inner], variables=(x, y)) for d in result.decompositions]
    assert all(
        (
            not field.contains(other.inner)
            for i, field in enumerate(fields)
            for j, other in enumerate(result.decompositions)
            if i != j
        )
    )


def test_factor_products():
    result = search_rational((x + y) ** 6, (x, y), parameter=t)
    assert result.candidate_count == 3 and len(result.decompositions) == 2
    assert len(result.factors) == 4
    assert any(d.degrees == (2, 3) for d in result.decompositions)
    assert_compositions(result)


def test_candidate_count():
    result = search_rational((x + y) ** 12, (x, y), parameter=t, max_candidates=1)
    factors, expected = (result.factors, 0)
    targets = {2, 3, 4, 6}
    comparison = result.comparison_variables
    for length in range(1, len(factors) + 1):
        for subset in combinations(factors, length):
            product = sp.Poly(sp.prod(p.as_expr() for p in subset), *comparison, x, y)
            dx = max(sum(m[:2]) for m in product.monoms())
            dy = max(sum(m[2:]) for m in product.monoms())
            diagonal = sp.Poly(
                product.as_expr().xreplace(dict(zip(comparison, (x, y), strict=True))), x, y
            )
            expected += int(dx == dy and dx in targets and diagonal.is_zero)
    assert result.candidate_count == expected and expected > 1
    assert result.examined == 1 and (not result.complete) and result.verify()


@pytest.mark.parametrize("function", [x**4 + y, (x**4 + y) / (x + y + 1), x**4 + y**4])
def test_negative_search(function):
    result = search_rational(function, (x, y), parameter=t)
    assert result.complete and result.status == "indecomposable" and (not result.decompositions)
    assert result.factors and result.verify()


@pytest.mark.parametrize("function", [0, 2, x + y, (x + y) / (x - y)])
def test_low_degree(function):
    result = search_rational(function, (x, y), parameter=t)
    assert result.status == "degree_below_two" and result.complete and result.verify()
    assert not result.factors and (not result.decompositions)


@pytest.mark.parametrize("constant", [sp.sqrt(2), sp.I])
def test_number_fields(constant):
    constants = sp.QQ.algebraic_field(constant)
    h = x * x + constant * x * y + y * y
    result = search_rational(h * h + constant, (x, y), constants=constants, parameter=t)
    assert result.complete and result.status == "decomposable"
    assert result.decompositions[0].constants == constants
    assert_compositions(result)


def test_visible_poles():
    function = (x * x + y * y) ** 2
    original = sp.Mul(function * (x - y), sp.Pow(x - y, -1, evaluate=False), evaluate=False)
    result = search_rational(original, (x, y), parameter=t)
    decomposition = result.decompositions[0]
    assert decomposition.original_expression == original
    assert sp.expand(decomposition.function - function) == 0
    assert decomposition.input_guard == sp.Ne(x - y, 0)
    assert_compositions(result)


@pytest.mark.parametrize("variables", [(x, y), (y, x), (z, y, x)])
def test_coordinate_order(variables):
    result = search_rational((x * x + y * y) ** 2, list(variables), parameter=t)
    assert result.variables == variables and result.parameter is t
    assert_compositions(result)
    a = sp.Dummy("x", positive=True)
    b = sp.Symbol("x", real=True)
    result = search_rational((a * a + b * b) ** 2, (a, b), parameter=t)
    assert result.variables[0] is a and result.variables[1] is b
    assert_compositions(result)


def test_iterable_inputs():
    result = search_rational(sp.Poly((x * x + y * y) ** 2, x, y), (v for v in (x, y)), parameter=t)
    assert_compositions(result)


@pytest.mark.parametrize("expression", [x + z, x + sp.sqrt(2), sp.sin(x), x + 0.5, sp.oo, sp.zoo])
def test_invalid_input(expression):
    with pytest.raises(DecompositionError):
        search_rational(expression, (x, y))


@pytest.mark.parametrize(
    "options",
    [
        {"variables": ()},
        {"variables": (x, x)},
        {"variables": (x, 1)},
        {"variables": (x, sp.Symbol("n", commutative=False))},
        {"parameter": x},
        {"parameter": 1},
        {"max_candidates": -1},
        {"max_candidates": True},
        {"max_degree": 0},
        {"max_degree": False},
        {"constants": sp.RR},
        {"constants": sp.GF(7)},
        {"constants": sp.QQ.frac_field(z)},
    ],
)
def test_invalid_options(options):
    kwargs = {"variables": (x, y)}
    kwargs.update(options)
    with pytest.raises(ValueError):
        search_rational((x * x + y * y) ** 2, **kwargs)


def test_altered_certificates():
    result = search_rational((x * x + y * y) ** 2, (x, y), parameter=t)
    decomposition = result.decompositions[0]
    for changes in (
        {"outer": t * t + 1},
        {"inner": x + y},
        {"outer": t * t + x},
        {"inner": x * x + z * z},
        {"parameter": x},
        {"variables": (x, x)},
        {"constants": sp.ZZ},
    ):
        assert not replace(decomposition, **changes).verify()
    for changes in (
        {"decompositions": ()},
        {"factors": result.factors[:-1]},
        {"factors": (result.factors[0] * 2, *result.factors[1:])},
        {"examined": -1},
        {"examined": True},
        {"examined": 0},
        {"candidate_count": result.candidate_count + 1},
        {"comparison_variables": (x, z)},
        {"comparison_variables": ()},
        {"comparison_variables": (z, z)},
        {"parameter": x},
        {"constants": sp.ZZ},
    ):
        assert not replace(result, **changes).verify()
    assert not RationalDecomposition(x + y, (x, y), t, sp.QQ, t, x + y).verify()


@pytest.mark.parametrize("function", [(x + y) ** 2, (x + y) ** 3])
def test_prime_degree(function):
    result = search_rational(function, (x, y))
    assert result.complete and result.status == "indecomposable"
    assert result.factors == () and result.candidate_count == 0
    assert result.decompositions == () and result.verify()


def test_degree_limit():
    from alggeom.decomposition import DegreeLimitError

    with pytest.raises(DegreeLimitError, match="total degree"):
        search_rational((x + y) ** 12, (x, y), max_degree=8)
