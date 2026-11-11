import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.modules import ModuleHomomorphism, PresentedModule, direct_sum, tensor_product

x, y, t = sp.symbols("x y t")


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(-3)), sp.QQ.frac_field(t)])
def test_module_maps(field):
    ring = PolynomialRing((x, y), field)
    source = PresentedModule(ring, [[x * x]])
    target = PresentedModule(ring, [[x]])
    projection = source.hom([[1]], target)
    assert projection.is_surjective and not projection.is_injective
    assert projection.cokernel.is_zero
    assert projection.kernel.annihilator() == ring.ideal(x)
    generator = projection.kernel_map(projection.kernel.gens[0])
    assert generator == source([x])
    assert projection(generator) == target.zero
    assert source([x * x]) == source.zero
    assert (source([1]) + source([x])) * x == source([x])
    with pytest.raises(TypeError):
        source([1]) + target([1])
    with pytest.raises(ValueError, match="relations"):
        target.hom([[1]], source)


def test_hom_and_tensor():
    ring = PolynomialRing((x, y))
    left = PresentedModule(ring, [[x]])
    right = PresentedModule(ring, [[x * x]])
    hom = left.hom_module(right)
    assert hom.annihilator() == ring.ideal(x)
    mapping = left.hom([[x]], right)
    element = hom.from_homomorphism(mapping)
    assert hom.homomorphism(element)(left.gens[0]) == right([x])
    assert hom.from_homomorphism(hom.homomorphism(element)) == element
    assert tensor_product(left, right).annihilator() == ring.ideal(x)
    other = PresentedModule(ring, [[y]])
    assert tensor_product(left, other).annihilator() == ring.ideal(x, y)
    assert direct_sum(left, other).annihilator() == ring.ideal(x * y)


def test_quotient_and_inverse():
    ring = PolynomialRing((x,)).quotient(x * x)
    module = PresentedModule(ring, sp.zeros(1, 0))
    action = module.hom([[x]], module)
    assert action.kernel_map.matrix == sp.Matrix([x])
    assert action.kernel.annihilator() == ring.ambient.ideal(x)
    assert action.image.annihilator() == ring.ambient.ideal(x)
    assert action.cokernel.annihilator() == ring.ambient.ideal(x)
    unit = module.hom([[1 + x]], module)
    inverse = unit.inverse()
    assert unit.compose(inverse)(module.gens[0]) == module.gens[0]
    assert inverse.matrix == sp.Matrix([1 - x])
    with pytest.raises(ValueError, match="isomorphism"):
        action.inverse()


def test_empty_modules():
    ring = PolynomialRing((x,))
    zero = PresentedModule(ring, sp.zeros(0, 0))
    free = PresentedModule(ring, sp.zeros(2, 0))
    assert tensor_product(zero, free).is_zero
    assert free.hom_module(zero).is_zero
    assert zero.hom_module(free).is_zero
    assert direct_sum(zero, free).generator_count == 2
    mapping = ModuleHomomorphism(zero, free, sp.zeros(2, 0))
    assert mapping.kernel.is_zero and mapping.image.is_zero
    assert not mapping.is_surjective
