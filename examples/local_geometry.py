"""Distinguish a cusp's local dimension, tangent dimension and intersection multiplicity."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.local import intersection_multiplicity
from alggeom.schemes import AffineScheme

x, y = sp.symbols("x y")
S = PolynomialRing((x, y))
cusp = AffineScheme(S.ideal(y**2 - x**3))
assert cusp.local_dimension((0, 0)) == 1
assert cusp.embedding_dimension((0, 0)) == 2
assert cusp.tangent_cone((0, 0)) == S.ideal(y**2)
assert cusp.multiplicity((0, 0)) == 2
assert intersection_multiplicity(cusp, S.ideal(y), (0, 0)) == 3
assert cusp.is_smooth((1, 1)) and not cusp.is_smooth((0, 0))
print(
    "Cusp: local dimension 1, tangent dimension 2, multiplicity 2; y=0 meets it with multiplicity 3."
)
