"""Bounded generated identities checked independently of certificate replay."""

import sympy as sp
from algroots.function_fields import RationalCoordinateMap
from hypothesis import given, settings
from hypothesis import strategies as st

from alggeom.decomposition.rational import search_rational

x, y, t = sp.symbols("x y t")


@given(st.integers(-3, 3), st.integers(-3, 3), st.sampled_from([(2, 2), (2, 3), (3, 2)]))
@settings(max_examples=18, deadline=None, derandomize=True)
def test_generated_compositions(a, b, degrees):
    outer_degree, inner_degree = degrees
    inner = x**inner_degree + a * x + b
    function = sp.expand(inner**outer_degree + inner)
    result = search_rational(function, x, parameter=t)
    assert result.complete
    answers = [d for d in result.decompositions if d.degrees == degrees]
    assert len(answers) == 1
    answer = answers[0]
    assert sp.cancel(answer.outer.subs(t, answer.inner) - function) == 0
    coordinates = RationalCoordinateMap(sp.QQ.frac_field(x), answer.inner, t)
    expressed = coordinates.express(inner)
    assert sp.cancel(expressed.subs(t, answer.inner) - inner) == 0
    reverse = RationalCoordinateMap(sp.QQ.frac_field(x), inner, t).express(answer.inner)
    assert sp.cancel(reverse.subs(t, inner) - answer.inner) == 0


@given(st.integers(-3, 3), st.integers(1, 4))
@settings(max_examples=10, deadline=None, derandomize=True)
def test_mobius_coordinates(shift, scale):
    inner = x * x + y * y
    transformed = sp.cancel((scale * inner + shift) / (inner + shift + scale + 1))
    field = RationalCoordinateMap(sp.QQ.frac_field(x, y), inner, t)
    expression = field.express(transformed)
    assert sp.cancel(expression.subs(t, inner) - transformed) == 0
    inverse = RationalCoordinateMap(sp.QQ.frac_field(x, y), transformed, t).express(inner)
    assert sp.cancel(inverse.subs(t, transformed) - inner) == 0


@given(st.sampled_from([(1, 1, 0, 1), (0, 1, 1, 0), (2, 1, 1, 1)]))
@settings(max_examples=3, deadline=None, derandomize=True)
def test_mobius_search(coefficients):
    from algroots.function_fields import luroth_generator

    a, b, c, d = coefficients
    original = (x + 1 / x) ** 2
    transformed = sp.cancel((a * original + b) / (c * original + d))
    result = search_rational(transformed, x, parameter=t)
    assert result.complete and len(result.decompositions) == 3
    expected = (x * x, x + 1 / x, x - 1 / x)
    for generator in expected:
        field = luroth_generator([generator], variables=(x,), parameter=t)
        matches = [pair for pair in result.decompositions if field.contains(pair.inner)]
        assert len(matches) == 1
        pair = matches[0]
        assert sp.cancel(pair.outer.subs(t, pair.inner) - transformed) == 0
        recovered = luroth_generator([pair.inner], variables=(x,), parameter=t)
        coordinates = recovered.express(generator)
        assert sp.cancel(coordinates.expression.subs(t, recovered.generator) - generator) == 0
