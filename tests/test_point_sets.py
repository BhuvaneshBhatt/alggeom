"""Independent complete point sets and exact candidate-budget boundaries."""

from itertools import product

import pytest
import sympy as sp
from algroots.finite_fields import field_elements

from alggeom.fields import field_order, finite_field
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y = sp.symbols("x y")


@pytest.mark.parametrize("field", [sp.GF(2), sp.GF(3), finite_field(2, 2)])
def test_complete_point_sets(field):
    q = field_order(field)
    pool = tuple(field_elements(field))
    C = AffineScheme(PolynomialRing((x, y), field).ideal(y - x * x))
    expected = {
        (field.to_sympy(a), field.to_sympy(b)) for a, b in product(pool, repeat=2) if b == a * a
    }
    assert {p.coordinates for p in C.rational_points(max_candidates=q * q)} == expected
    with pytest.raises(RuntimeError, match="max_candidates"):
        C.rational_points(max_candidates=q * q - 1)
    P = ProjectiveScheme(PolynomialRing((x, y), field))
    expected = {(field.to_sympy(field.one), field.to_sympy(a)) for a in pool} | {
        (sp.S.Zero, sp.S.One)
    }
    assert {p.coordinates for p in P.rational_points(max_candidates=q + 1)} == expected
    with pytest.raises(RuntimeError, match="max_candidates"):
        P.rational_points(max_candidates=q)
