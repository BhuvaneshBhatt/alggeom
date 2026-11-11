from dataclasses import replace

import pytest
import sympy as sp

from alggeom.fields import finite_field
from alggeom.ideals import PolynomialRing
from alggeom.noether import noether_normalization

x, y, z = sp.symbols("x y z")


@pytest.mark.parametrize(
    "field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(-3)), sp.GF(2), finite_field(3, 2)]
)
def test_finite_projection(field):
    ideal = PolynomialRing((x, y), field).ideal(x * y - 1)
    result = noether_normalization(ideal)
    assert result.verify()
    assert len(result.base.variables) == 1 and result.generic_rank == 2
    assert result.base_embedding.is_injective
    assert len(result.integral_relations) == 1
    assert result.projection.source.coordinate_ring == result.ring
    assert (
        result.inverse_map.compose(result.coordinate_map).images
        == result.coordinate_map.source.gens
    )
    assert not replace(result, integral_relations=(sp.S.Zero,)).verify()
    assert not replace(result, module_generators=()).verify()
    assert not replace(result, module_generators=result.module_generators[:-1]).verify()


def test_finite_field_nonlinear_change():
    ideal = PolynomialRing((x, y), sp.GF(2)).ideal(x * x * y + x * y * y + 1)
    result = noether_normalization(ideal, max_linear_attempts=2)
    assert result.verify() and result.generic_rank > 1
    with pytest.raises(RuntimeError, match="max_transform_degree"):
        noether_normalization(ideal, max_linear_attempts=2, max_transform_degree=2)


def test_module_spanning():
    ideal = PolynomialRing((x, y, z)).ideal(x * y, x * z)
    result = noether_normalization(ideal)
    assert len(result.base.variables) == 2
    assert result.generic_rank == 1
    assert len(result.module_generators) > result.generic_rank
    assert result.verify()
    artinian = noether_normalization(PolynomialRing((x, y)).ideal(x * x, y * y))
    assert artinian.base.variables == () and artinian.generic_rank == 4
    with pytest.raises(ValueError, match="nonzero"):
        noether_normalization(PolynomialRing((x,)).ideal(1))


def test_projection_binding():
    result = noether_normalization(PolynomialRing((x, y)).ideal(x * y - 1))
    parameter = result.base.variables[0]
    squared = result.base.hom((parameter**2,), result.ring)
    assert squared.is_injective
    generic_value = sp.Dummy("generic_value")
    generic_fiber = PolynomialRing(result.ring.variables, sp.QQ.frac_field(generic_value)).quotient(
        *result.ring.ideal.generators, parameter**2 - generic_value
    )
    assert generic_fiber.vector_space_dimension == 4
    assert not replace(result, base_embedding=squared).verify()
    unrelated = result.base.hom((parameter,), result.base)
    assert not replace(result, base_embedding=unrelated).verify()
    source = result.coordinate_map.source
    assert not replace(result, inverse_map=source.hom(source.variables, source)).verify()


@pytest.mark.parametrize("field", [sp.QQ, sp.GF(2), finite_field(3, 2)])
def test_generic_rank_lex_localization(field, monkeypatch):
    from algroots.quotient import QuotientAlgebra

    ideal = PolynomialRing((x, y, z), field).ideal(x * y, x * z)
    normalization = noether_normalization(ideal)
    count = len(normalization.ring.variables) - len(normalization.base.variables)
    fraction_field = field.frac_field(*normalization.base.variables)
    reference = QuotientAlgebra.from_polynomials(
        normalization.ring.ideal.generators,
        normalization.ring.variables[:count],
        domain=fraction_field,
        groebner_backend="sympy",
    )
    monkeypatch.setattr(
        QuotientAlgebra, "from_polynomials", lambda *a, **k: pytest.fail("recomputed generic basis")
    )
    assert normalization.generic_rank == reference.dimension == 1
