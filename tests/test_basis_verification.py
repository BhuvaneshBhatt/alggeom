"""Independent basis-certificate checks for candidates supplied by algroots."""

import pytest
import sympy as sp
from algroots.groebner import compute_groebner_basis
from algroots.native_f4 import compiled_available
from sympy.polys.rings import PolyRing

from alggeom.groebner_certificates import verify_reduced_groebner_basis
from alggeom.ideals import PolynomialRing


@pytest.mark.parametrize("order", ("lex", "grlex", "grevlex"))
def test_term_orders(order):
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z), sp.QQ, order)
    ideal = ring.ideal(x * y - z, x**2 - y, y**2 - z)
    sparse_ring = PolyRing((x, y, z), sp.QQ, order)
    candidate = ideal.groebner_basis_from("sympy")
    assert candidate == ideal.groebner_basis
    assert verify_reduced_groebner_basis(ideal.generators, candidate, sparse_ring)
    assert not verify_reduced_groebner_basis(ideal.generators, candidate[:-1], sparse_ring)


def test_candidate_tampering():
    x, y = sp.symbols("x y")
    sparse_ring = PolyRing((x, y), sp.QQ, "lex")
    assert not verify_reduced_groebner_basis((x, y), (x,), sparse_ring)
    assert not verify_reduced_groebner_basis((x, y), (x + y, y + 1), sparse_ring)
    # Generates the same ideal but is not reduced.
    assert not verify_reduced_groebner_basis((x, y), (x + y, y), sparse_ring)
    assert verify_reduced_groebner_basis((), (), sparse_ring)
    assert not verify_reduced_groebner_basis((x,), (), sparse_ring)
    with pytest.raises(ValueError, match="coefficient field"):
        verify_reduced_groebner_basis((x,), (x,), PolyRing((x, y), sp.ZZ, "lex"))


def _algroots_candidate(expressions, variables, order):
    return compute_groebner_basis(
        expressions, PolyRing(variables, sp.QQ, order), backend="native", certify=False
    )


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_algroots_candidate():
    x, y, z = sp.symbols("x y z")
    equations = (x * y - z, x**2 - y, y**2 - z)
    for order in ("grevlex", "lex"):
        ring = PolyRing((x, y, z), sp.QQ, order)
        generators = tuple(ring.from_expr(expression) for expression in equations)
        candidate = _algroots_candidate(equations, (x, y, z), order)
        assert verify_reduced_groebner_basis(generators, candidate, ring)


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_rational_candidate():
    x, y = sp.symbols("x y")
    equations = (x**2 / 2 + y / 3, x * y / 5 - 1)
    ring = PolyRing((x, y), sp.QQ, "lex")
    generators = tuple(ring.from_expr(expression) for expression in equations)
    candidate = _algroots_candidate(equations, (x, y), "lex")
    assert verify_reduced_groebner_basis(generators, candidate, ring)
