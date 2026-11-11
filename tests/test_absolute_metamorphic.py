"""Metamorphic tests of exact primary decomposition under field and coordinate changes."""

import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing


def test_primary_coordinate_transport():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    original = absolute_primary_decomposition(ring.ideal(x**2, x * y))
    transformed = absolute_primary_decomposition(ring.ideal((x + y) ** 2, (x + y) * y))
    expected = tuple(
        ring.ideal(*(g.subs(x, x + y) for g in primary.generators))
        for primary in original.components
    )
    assert len(transformed.components) == len(expected)
    assert all(component in transformed.components for component in expected)
    assert sorted(c.krull_dimension for c in transformed.components) == sorted(
        c.krull_dimension for c in original.components
    )


def test_nilpotents_after_field_extension():
    x, y = sp.symbols("x y")
    rational = PolynomialRing((x, y))
    extension = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.I))
    original = absolute_primary_decomposition(rational.ideal(x**2, x * y))
    extended = absolute_primary_decomposition(extension.ideal(x**2, x * y))
    assert len(original.components) == len(extended.components) == 2
    assert sorted(c.krull_dimension for c in original.components) == sorted(
        c.krull_dimension for c in extended.components
    )
    for primary in original.components:
        assert primary.change_ring(extension) in extended.components


def test_graded_coordinate_transport():
    from alggeom.graded import degree, hilbert_function, hilbert_series
    from alggeom.local import multiplicity

    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    ideal = ring.ideal(x**2, x * y)
    changed = ring.ideal((x + y) ** 2, (x + y) * y)
    assert degree(ideal) == degree(changed) == 1
    assert multiplicity(ideal, (0, 0)) == multiplicity(changed, (0, 0)) == 1
    assert sp.cancel(hilbert_series(ideal) - hilbert_series(changed)) == 0
    for n in range(7):
        assert hilbert_function(ideal, n) == hilbert_function(changed, n)
    original = absolute_primary_decomposition(ideal)
    transformed = absolute_primary_decomposition(changed)
    for component in original.components:
        corresponding = ring.ideal(*(g.subs(x, x + y) for g in component.generators))
        assert corresponding in transformed.components
        assert degree(component) == degree(corresponding)


def test_term_order_transport():
    x, y = sp.symbols("x y")
    before = PolynomialRing((x, y), order="lex").ideal(x * x + y, x * y)
    after = PolynomialRing((y, x), order="lex").ideal(y * y + x, y * x)
    swap = {x: y, y: x}
    expected = tuple(sp.expand(g.subs(swap, simultaneous=True)) for g in before.groebner_basis)
    # Transport the variable ordering as well as the ideal. A general linear
    # coordinate change would not preserve the initial monomial ideal.
    assert all(after.contains(g) for g in expected)
    assert all(before.contains(g.subs(swap, simultaneous=True)) for g in after.groebner_basis)
    expected_initial = sorted(
        sp.Poly(g, x, y).LM(order="lex").exponents for g in before.groebner_basis
    )
    actual_initial = sorted(
        sp.Poly(g, y, x).LM(order="lex").exponents for g in after.groebner_basis
    )
    assert expected_initial == actual_initial


def test_quadratic_splitting_under_field_extension():
    x, y = sp.symbols("x y")
    rational = PolynomialRing((x, y))
    complex_field = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.I))
    source = rational.ideal(x**2 + y**2)
    split = absolute_primary_decomposition(source)
    assert len(split.components) == 2
    factors = (complex_field.ideal(x + sp.I * y), complex_field.ideal(x - sp.I * y))
    assert all(
        any(part == factor.change_ring(part.ring) for part in split.components)
        for factor in factors
    )
    assert factors[0].intersection(factors[1]) == source.change_ring(complex_field)
    assert all(part.krull_dimension == 1 for part in factors)
