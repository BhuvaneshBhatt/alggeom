"""Certified nonprincipal decomposition and recursive exceptional-locus tests."""

import pytest
import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.groebner_cover import comprehensive_groebner_cover
from alggeom.ideals import PolynomialRing


def test_disjoint_quadratic_components():
    x, y, u, v = sp.symbols("x y u v")
    ring = PolynomialRing((x, y, u, v))
    ideal = ring.ideal(x**2 + y**2, u**2 + v**2)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 4
    extended = ideal.change_ring(decomposition.components[0].ring)
    reconstructed = extended.ring.ideal(1)
    for component in decomposition.components:
        assert component.krull_dimension == 2
        assert len(component.groebner_basis) == 2
        reconstructed = reconstructed.intersection(component)
    assert reconstructed == extended


def test_disjoint_nonreduced_primary_components():
    x, y, u, v = sp.symbols("x y u v")
    ring = PolynomialRing((x, y, u, v))
    ideal = ring.ideal((x**2 + y**2) ** 2, u**2 + v**2)
    result = absolute_primary_decomposition(ideal)
    assert len(result.components) == 4
    product = result.components[0].ring.ideal(1)
    for component in result.components:
        product = product.intersection(component)
    assert product == ideal.change_ring(product.ring)
    assert all(
        any(sp.Poly(g, x, y, u, v).total_degree() > 1 for g in component.generators)
        for component in result.components
    )


def test_nested_nonlin_exceptional_cover():
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover(((a**2 + b**2) * x, a * y), (x, y), (a, b))
    expected = {
        (0, 0): (),
        (0, 1): (x,),
        (1, 0): (x, y),
        (sp.I, 1): (y,),
        (-sp.I, 1): (y,),
    }
    for values, basis in expected.items():
        assert cover.verify_at(values)
        assert tuple(sp.expand(item) for item in cover.specialize(values)) == basis


def test_intersection_of_nonlinear_graphs():
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover(((a**2 + b**2) * x, b * y), (x, y), (a, b))
    for values in ((0, 0), (0, 1), (1, 0), (sp.I, 1), (-sp.I, 1), (2, 3)):
        assert cover.verify_at(values)


def test_geometric_integrality_disjoint_blocks():
    from alggeom.geometric import geometric_components, is_geometrically_integral

    x, y, u, v = sp.symbols("x y u v")
    ring = PolynomialRing((x, y, u, v))
    irreducible = ring.ideal(x**2 + y**2 + 1, u**2 + v**2 + 1)
    assert is_geometrically_integral(irreducible)
    split = ring.ideal(x**2 + y**2, u**2 + v**2)
    assert not is_geometrically_integral(split)
    components = geometric_components(split)
    assert len(components.components) == 4


def test_component_budget():
    from alggeom.components import DecompositionLimitError

    x, y, u, v = sp.symbols("x y u v")
    ideal = PolynomialRing((x, y, u, v)).ideal(x**2 + y**2, u**2 + v**2)
    with pytest.raises(DecompositionLimitError, match="max_components"):
        absolute_primary_decomposition(ideal, max_components=3)
    with pytest.raises(ValueError, match="max_components"):
        absolute_primary_decomposition(ideal, max_components=0)


def test_cover_piece_budget():
    from alggeom.groebner_cover import GroebnerCoverLimitError

    a, b, x, y = sp.symbols("a b x y")
    with pytest.raises(GroebnerCoverLimitError, match="piece budget"):
        comprehensive_groebner_cover(((a**2 + b**2) * x, a * y), (x, y), (a, b), max_pieces=2)


def test_polynomial_graph_geometric_prime():
    from alggeom.geometric import geometric_components, is_geometrically_integral

    x, y, z, w = sp.symbols("x y z w")
    ideal = PolynomialRing((x, y, z, w)).ideal(z - x * y, w - x**2)
    decomposition = absolute_primary_decomposition(ideal)
    assert decomposition.components == (ideal,)
    assert is_geometrically_integral(ideal)
    assert geometric_components(ideal).components == (ideal,)


def test_graph_criterion_rejects_extra_relation():
    from alggeom.geometric import _graph_prime_decomposition

    x, y, z, w = sp.symbols("x y z w")
    ideal = PolynomialRing((x, y, z, w)).ideal(z - x * y, w - x**2, x * y)
    assert _graph_prime_decomposition(ideal) is None
