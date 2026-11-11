"""Primary components under certified polynomial graph equivalences."""

import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing


def test_graph_hypersurface_splitting():
    x, y, z, w = sp.symbols("x y z w")
    ring = PolynomialRing((x, y, z, w))
    ideal = ring.ideal(x**2 + y**2, z - x * y)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 2
    recovered = decomposition.components[0].ring.ideal(1)
    for component in decomposition.components:
        assert component.krull_dimension == 2
        recovered = recovered.intersection(component)
    assert recovered == ideal.change_ring(recovered.ring)
    assert all(component.contains(z - x * y) for component in decomposition.components)
    assert {component.contains(x + sp.I * y) for component in decomposition.components} == {
        True,
        False,
    }


def test_graph_embedded_primary_components():
    x, y, z, w = sp.symbols("x y z w")
    ideal = PolynomialRing((x, y, z, w)).ideal(x**2 + 1, y**2, y * z, w - x * y)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 4
    assert sorted(component.krull_dimension for component in decomposition.components) == [
        0,
        0,
        1,
        1,
    ]
    recovered = decomposition.components[0].ring.ideal(1)
    for component in decomposition.components:
        recovered = recovered.intersection(component)
    assert recovered == ideal.change_ring(recovered.ring)
