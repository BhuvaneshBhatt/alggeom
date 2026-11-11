"""Selected real graph equivalence for recursive composition chains."""

from dataclasses import replace

import pytest
import sympy as sp
from semialg.function_verification import RealCompositionVerification, RealFunctionVerification

from alggeom.decomposition.chains import algebraic_chain, rational_chain
from alggeom.decomposition.rational import RationalDecomposition
from alggeom.decomposition.real import verify_chain

x, y, z, u, v = sp.symbols("x y z u v", real=True)


def test_rational_graph():
    chain = rational_chain(x**8, x)
    result = verify_chain(chain, output=z, intermediates=(u, v))
    assert type(result) is RealCompositionVerification
    assert all(type(check) is RealFunctionVerification for check in result.function_checks)
    assert result.equivalent and result.certified and result.function_verified
    assert result.comparison.variables == (x, z)
    assert result.composition_formula.free_symbols == {x, z}


def test_wrong_branch():
    chain = rational_chain((x * x + y * y) ** 4, (x, y))
    result = verify_chain(
        chain,
        output=z,
        intermediates=(u, v),
        component_selectors=(sp.true, sp.true, v < 0),
        check_functions=False,
    )
    assert result.composition_formula is sp.false
    assert not result.equivalent and result.certified and (result.function_verified is None)
    assert result.comparison.variables == (x, y, z)


def test_incomplete_chain():
    chain = rational_chain(x**8, x, max_depth=1)
    assert not chain.complete
    result = verify_chain(chain, output=z)
    assert result.equivalent and result.certified and result.function_verified


def test_original_hole():
    original = sp.Mul(x**4 * (x - 1), sp.Pow(x - 1, -1, evaluate=False), evaluate=False)
    chain = rational_chain(original, x)
    result = verify_chain(chain, output=z, check_functions=False)
    assert not result.equivalent and result.certified
    restricted = verify_chain(chain, output=z, domain=sp.Ne(x, 1))
    assert restricted.equivalent and restricted.certified and restricted.function_verified


def test_principal_root():
    chain = algebraic_chain(z**8 - x, (z, x))
    assert len(chain.components) == 3 and chain.verify()
    result = verify_chain(
        chain, domain=x >= 0, target_selector=z >= 0, component_selectors=(z >= 0, sp.true, sp.true)
    )
    assert result.equivalent and result.certified and result.function_verified
    wrong = verify_chain(
        chain,
        domain=x > 0,
        target_selector=z >= 0,
        component_selectors=(z < 0, sp.true, sp.true),
        check_functions=False,
    )
    assert not wrong.equivalent and wrong.certified


def test_internal_pole():
    initial = RationalDecomposition(1 / (1 + x**8), (x,), u, sp.QQ, u**4 / (u**4 + 1), 1 / x**2)
    chain = rational_chain(initial.original_expression, x, initial=initial)
    assert chain.verify() and len(chain.components) == 3
    unrestricted = verify_chain(chain, output=z, check_functions=False)
    assert not unrestricted.equivalent and unrestricted.certified
    restricted = verify_chain(chain, output=z, domain=sp.Ne(x, 0))
    assert restricted.equivalent and restricted.certified and restricted.function_verified


def test_selector_sign():
    initial = RationalDecomposition(1 / x**6, (x,), u, sp.QQ, u**2, 1 / x**3)
    chain = rational_chain(initial.original_expression, x, initial=initial)
    positive = verify_chain(
        chain,
        output=z,
        intermediates=(u,),
        domain=x > 0,
        component_selectors=(sp.true, u > 0),
        check_functions=False,
    )
    negative = verify_chain(
        chain,
        output=z,
        intermediates=(u,),
        domain=x < 0,
        component_selectors=(sp.true, u > 0),
        check_functions=False,
    )
    assert positive.equivalent and positive.certified
    assert not negative.equivalent and negative.certified
    assert positive.eliminations[1].method == "guarded-rational-chain-substitution"


@pytest.mark.parametrize(
    "options",
    [
        {"output": x},
        {"intermediates": (u,)},
        {"intermediates": (u, u)},
        {"intermediates": (u, v), "component_selectors": (sp.true, sp.true, z > 0)},
        {"component_selectors": (sp.true,)},
        {"target_selector": u > 0},
        {"domain": z > 0},
    ],
)
def test_coordinate_scope(options):
    chain = rational_chain(x**8, x)
    kwargs = {"output": z, **options}
    with pytest.raises(ValueError):
        verify_chain(chain, **kwargs)


def test_constant_chains():
    for expression in [x + 1, sp.Integer(3)]:
        result = verify_chain(rational_chain(expression, x), output=z)
        assert result.equivalent and result.certified and result.function_verified


def test_invalid_chain():
    chain = rational_chain(x**4, x)
    bad = replace(chain, root=replace(chain.root, inner=None))
    with pytest.raises(ValueError, match="failed verification"):
        verify_chain(bad)
    K = sp.QQ.algebraic_field(sp.I)
    complex_chain = rational_chain(x**4 + sp.I, x, constants=K)
    with pytest.raises(ValueError, match="real embedding"):
        verify_chain(complex_chain)
    with pytest.raises(ValueError, match="chain is required"):
        verify_chain(chain.root.split)


def test_fixed_coordinates():
    chain = algebraic_chain(z - x**4, (z, x))
    with pytest.raises(ValueError, match="output variable"):
        verify_chain(chain, output=y)
    with pytest.raises(ValueError, match="intermediate variables"):
        verify_chain(chain, intermediates=(u,))
