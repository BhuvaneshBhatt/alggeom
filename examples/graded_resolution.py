"""The twisted cubic's three quadrics have two linear syzygies."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.resolutions import free_resolution

x, y, z, w, t = sp.symbols("x y z w t")
ring = PolynomialRing((x, y, z, w))
ideal = ring.ideal(x * z - y * y, x * w - y * z, y * w - z * z)
resolution = free_resolution(ideal)
assert resolution.shifts == ((0,), (2, 2, 2), (3, 3))
assert resolution.betti_table == {0: (1, 0, 0), 1: (0, 3, 2)}
assert resolution.verify()
assert sp.cancel(resolution.hilbert_series(parameter=t) - (1 + 2 * t) / (1 - t) ** 2) == 0
print("Betti rows (j-i), columns (i):", resolution.betti_table)
print(
    "Projective dimension:", resolution.projective_dimension, "regularity:", resolution.regularity
)
