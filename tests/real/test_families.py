"""Complete parameter strata and independent real function obligations."""

from dataclasses import replace

import pytest
import sympy as sp
from semialg.function_verification import verify_real_composition, verify_real_function

from alggeom.decomposition.parameters import polynomial_family

x, y, t, a, b, c = sp.symbols("x y t a b c", real=True)


def test_degree_strata():
    R = polynomial_family(a * x**4 + b * x * x, x, parameters=(a, b), component_variable=t)
    assert R.verify()
    cases = [
        ({a: 2, b: 3}, 4, "decomposable"),
        ({a: 0, b: 3}, 2, "indecomposable"),
        ({a: 0, b: 0}, -1, "degree_below_two"),
    ]
    for assignment, degree, status in cases:
        case = R.select(assignment)
        assert (case.degree, case.status) == (degree, status)
        assert case.polynomial == (a * x**4 + b * x * x).subs(assignment)
        for component in case.components:
            assert component.condition is sp.true
            assert sp.expand(component.composition - case.polynomial) == 0
    assert R.select({"a": 2, "b": 3}).components[0].outer == 2 * t * t + 3 * t
    with pytest.raises(ValueError):
        R.select({a: sp.Float(1), b: 2})


def test_undefined_strata():
    R = polynomial_family((a * x**4 + b * x * x) / c, x, parameters=(a, b, c))
    assert R.verify()
    assert R.select({a: 1, b: 1, c: 0}).status == "undefined"
    case = R.select({a: 1, b: 2, c: 3})
    assert case.status == "decomposable" and case.components[0].condition is sp.true
    assert sp.expand(case.components[0].composition - case.polynomial) == 0
    restricted = polynomial_family(x**4 / a, x, parameters=(a,), parameter_domain=a > 0)
    assert restricted.verify()
    assert restricted.select({a: 2}).status == "decomposable"
    with pytest.raises(ValueError, match="outside"):
        restricted.select({a: 0})


def test_exceptional_equation():
    R = polynomial_family(x**4 + a * x, x, parameters=(a,))
    assert R.verify()
    assert R.select({a: 0}).status == "decomposable"
    assert R.select({a: 1}).status == "indecomposable"
    yes = next(branch for branch in R.branches if branch.value.components)
    component = yes.value.components[0]
    bad_component = replace(component, outer=component.outer + 1)
    bad_branch = replace(yes, value=replace(yes.value, components=(bad_component,)))
    bad = replace(
        R,
        stratification=replace(
            R.stratification, branches=tuple(bad_branch if q is yes else q for q in R.branches)
        ),
    )
    assert not bad.verify()
    dropped = replace(R, stratification=replace(R.stratification, branches=(yes,)))
    assert not dropped.verify()


def test_ordered_degrees():
    R = polynomial_family(a * x**6 + b * x * x, x, parameters=(a, b))
    assert R.verify()
    assert {c.system.degrees for c in R.select({a: 1, b: 0}).components} == {(3, 2), (2, 3)}
    assert {c.system.degrees for c in R.select({a: 1, b: 1}).components} == {(3, 2)}
    assert R.select({a: 0, b: 1}).degree == 2


def test_disjunctive_domain():
    R = polynomial_family(
        a * x**4 + b * x * x, x, parameters=(a, b), parameter_domain=sp.Or(sp.Eq(a, 0), sp.Eq(b, 0))
    )
    assert R.verify()
    assert R.select({a: 2, b: 0}).degree == 4
    assert R.select({a: 0, b: 2}).degree == 2


@pytest.mark.parametrize(
    "f,degree", [(sp.Integer(0), -1), (sp.Integer(3), 0), (x + 1, 1), (x**3 + x, 3)]
)
def test_low_degree(f, degree):
    R = polynomial_family(f, x)
    assert R.verify() and R.select({}).degree == degree
    assert R.select({}).components == ()


def test_function_obligations():
    sqrt_graph = sp.And(sp.Eq(y * y, x), y >= 0)
    principal = verify_real_function(sqrt_graph, input=x, output=y, domain=x >= 0)
    assert principal.total and principal.single_valued and principal.certified
    partial = verify_real_function(sqrt_graph, input=x, output=y)
    assert not partial.total and partial.single_valued and partial.certified
    multivalued = verify_real_function(sp.Eq(y * y, x), input=x, output=y, domain=x > 0)
    assert multivalued.total and (not multivalued.single_valued)
    empty = verify_real_function(sp.false, input=x, output=y)
    assert not empty.total and empty.single_valued
    vacuous = verify_real_function(sp.false, input=x, output=y, domain=sp.false)
    assert vacuous.total and vacuous.single_valued


def test_parameter_totality():
    graph = sp.Eq(x * y, 1)
    unrestricted = verify_real_function(graph, input=x, output=y)
    assert not unrestricted.total and unrestricted.single_valued
    restricted = verify_real_function(graph, input=x, output=y, domain=sp.Ne(x, 0))
    assert restricted.total and restricted.single_valued
    parametric = verify_real_function(
        sp.Eq(a * y, x), input=x, output=y, parameters=(a,), domain=sp.Ne(a, 0)
    )
    assert parametric.total and parametric.single_valued


def test_multivalued_graphs():
    relation = sp.Or(sp.Eq(y, x), sp.Eq(y, -x))
    R = verify_real_composition(
        relation,
        sp.Or(sp.Eq(y, t), sp.Eq(y, -t)),
        sp.Eq(t, x),
        input=x,
        output=y,
        intermediate=t,
        check_functions=True,
    )
    assert R.equivalent and R.certified and (R.function_verified is False)
    S = verify_real_composition(
        sp.Eq(y, x + 1),
        sp.Eq(y, t + 1),
        sp.Eq(t, x),
        input=x,
        output=y,
        intermediate=t,
        check_functions=True,
    )
    assert S.equivalent and S.function_verified and S.certified


@pytest.mark.parametrize(
    "kwargs",
    [
        {"parameters": (a, a)},
        {"parameters": (x,)},
        {"parameter_domain": a > 0},
        {"component_variable": x},
    ],
)
def test_invalid_partition(kwargs):
    with pytest.raises(ValueError):
        polynomial_family(x**4, x, **kwargs)


@pytest.mark.parametrize("f", [1 / x, x**4 + sp.Float(0.5), x**4 + sp.I, x**4 + a])
def test_invalid_coefficients(f):
    with pytest.raises(ValueError):
        polynomial_family(f, x)


def test_graph_contract():
    with pytest.raises(ValueError, match="domain"):
        verify_real_function(sp.Eq(y, x), input=x, output=y, domain=y > 0)
    with pytest.raises(ValueError, match="undeclared"):
        verify_real_function(sp.Eq(y, x + a), input=x, output=y)


def test_graph_auxiliaries():
    from semialg.function_graph import semialgebraic_function_graph

    graph = semialgebraic_function_graph(sp.sqrt(x), y)
    result = verify_real_function(graph, input=x, output=y, domain=x >= 0)
    assert result.total and result.single_valued and result.certified
    assert not set(graph.auxiliary_variables) & result.formula.free_symbols


def test_empty_domain():
    domain = sp.And(sp.Eq(a * b - 1, 0), sp.Eq(b, 0))
    result = polynomial_family(x**4 / a, x, parameters=(a, b), parameter_domain=domain)
    assert result.verify()


def test_exact_parameter_names():
    family = polynomial_family(a * x**4 + b * x * x, x, parameters=(a, b))
    by_name = family.select({"a": 1, "b": 2})
    by_symbol = family.select({a: 1, b: 2})
    assert by_name == by_symbol
    assert by_name.polynomial == x**4 + 2 * x * x
    for assignment in ({"unknown": 1}, {a: sp.Float(1), b: 2}, {c: 1}):
        with pytest.raises(ValueError):
            family.select(assignment)
