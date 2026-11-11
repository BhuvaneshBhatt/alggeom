"""Composition of graph objects with bound auxiliaries and explicit domains."""

import pytest
import sympy as sp
from semialg.function_graph import semialgebraic_function_graph
from semialg.function_verification import verify_real_composition

x, y, z, t, a = sp.symbols("x y z t a", real=True)


def test_graph_auxiliaries():
    target = semialgebraic_function_graph(z ** sp.Rational(1, 4), x)
    outer = semialgebraic_function_graph(sp.sqrt(y), x)
    inner = semialgebraic_function_graph(sp.sqrt(z), y)
    result = verify_real_composition(
        target, outer, inner, input=z, output=x, intermediate=y, domain=z >= 0, check_functions=True
    )
    assert result.equivalent and result.certified and result.function_verified
    assert result.composition_formula.free_symbols == {x, z}


def test_graph_poles():
    target = semialgebraic_function_graph(z, x)
    outer = semialgebraic_function_graph(1 / y, x)
    inner = semialgebraic_function_graph(1 / z, y)
    unrestricted = verify_real_composition(target, outer, inner, input=z, output=x, intermediate=y)
    assert not unrestricted.equivalent
    restricted = verify_real_composition(
        target, outer, inner, input=z, output=x, intermediate=y, domain=sp.Ne(z, 0)
    )
    assert restricted.equivalent and restricted.certified


def test_graph_parameters():
    result = verify_real_composition(
        semialgebraic_function_graph(z + a, x),
        semialgebraic_function_graph(y + a, x),
        semialgebraic_function_graph(z, y),
        input=z,
        output=x,
        intermediate=y,
        parameters=(a,),
        domain=a >= 0,
    )
    assert result.equivalent and result.certified
    assert result.target_formula.free_symbols == {a, x, z}


@pytest.mark.parametrize(
    "graphs,options,reason",
    [
        ((sp.Eq(x, z), sp.Eq(x, y + z), sp.Eq(y, z)), {}, "undeclared"),
        ((sp.Eq(x, z), sp.Eq(x, y), sp.Eq(y, z)), {"domain": x > 0}, "domain"),
        ((semialgebraic_function_graph(z, t), sp.Eq(x, y), sp.Eq(y, z)), {}, "target"),
    ],
)
def test_graph_scope(graphs, options, reason):
    with pytest.raises(ValueError, match=reason):
        verify_real_composition(*graphs, input=z, output=x, intermediate=y, **options)
