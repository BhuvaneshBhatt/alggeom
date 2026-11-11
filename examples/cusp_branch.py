"""Normalize a cusp, retain its finite fiber, and select the upper real branch."""

import sympy as sp
from semialg.function_verification import verify_real_function

from alggeom.ideals import PolynomialRing
from alggeom.normalization import integral_closure
from alggeom.schemes import AffineScheme

x, y, t = sp.symbols("x y t", real=True)
S = PolynomialRing((x, y))
cusp = AffineScheme(S.ideal(y * y - x**3))
closure = integral_closure(cusp)
assert closure.verify() and closure.mapping.is_injective
K = cusp.coordinate_ring.fraction_field()
parameter = K(y / x)
assert parameter**2 == K(x) and parameter**3 == K(y)
fiber = cusp.coordinate_ring.quotient(x - 1).finite_algebra()
assert fiber.dimension == fiber.geometric_solution_count == 2
assert {p.coordinates for p in AffineScheme(S.ideal(y * y - x**3, x - 1)).rational_points()} == {
    (1, -1),
    (1, 1),
}
selected = sp.And(sp.Eq(y * y, x**3), x >= 0, y >= 0)

branch = verify_real_function(selected, input=x, output=y, domain=x >= 0)
assert branch.total and branch.single_valued and branch.certified
assert bool(selected.subs({x: 1, y: 1})) and not bool(selected.subs({x: 1, y: -1}))
# On x>0 the rational inverse y/x already proves both coordinate identities;
# the missing point has the explicit preimage t=0.
assert (t * t).subs(t, 0) == (t**3).subs(t, 0) == 0
print("The normalization has parameter y/x; x=1 has two points and y>=0 selects (1,1).")
