"""Exact sequences and tensor/Hom evaluation, independent of chosen syzygy matrices."""

import pytest
import sympy as sp

from alggeom.fields import finite_field
from alggeom.ideals import PolynomialRing
from alggeom.modules import PresentedModule, direct_sum, tensor_product

x, y = sp.symbols("x y")


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(2)), finite_field(3, 2)])
def test_short_exact_sequence(field):
    S = PolynomialRing((x, y), field)
    M = PresentedModule(S, [[x * x]])
    N = PresentedModule(S, [[x]])
    quotient = M.hom([[1]], N)
    inclusion = N.hom([[x]], M)
    assert inclusion.is_injective and quotient.is_surjective
    assert quotient.compose(inclusion)(N.gens[0]) == N.zero
    assert quotient.kernel.annihilator() == N.annihilator()
    for value in (0, 1, y, 1 + y):
        element = inclusion(N([value]))
        assert quotient(element) == N.zero
        assert quotient.kernel_map(quotient.kernel_map.lift(element)) == element


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(2)), finite_field(3, 2)])
def test_direct_sum_maps(field):
    S = PolynomialRing((x, y), field)
    M = PresentedModule(S, [[x]])
    N = PresentedModule(S, [[y]])
    total = direct_sum(M, N)
    i = M.hom([[1], [0]], total)
    j = N.hom([[0], [1]], total)
    p = total.hom([[1, 0]], M)
    q = total.hom([[0, 1]], N)
    for vector in ([1, 0], [0, 1], [x + y, 1 + y]):
        element = total(vector)
        assert i(p(element)) + j(q(element)) == element


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(2)), finite_field(3, 2)])
def test_tensor_and_hom(field):
    S = PolynomialRing((x,), field)
    M = PresentedModule(S, [[x]])
    N = PresentedModule(S, [[x * x]])
    product = tensor_product(M, N)
    assert product([x]) == product.zero
    H = M.hom_module(N)
    f = M.hom([[x]], N)
    h = H.from_homomorphism(f)
    for a, b in ((1, 2), (x + 1, 3), (2, 1 - x)):
        assert H.homomorphism((a + b) * h)(M.gens[0]) == a * f(M.gens[0]) + b * f(M.gens[0])
