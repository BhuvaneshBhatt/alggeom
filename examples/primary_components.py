"""An embedded component is invisible in the reduced support."""

from functools import reduce

import sympy as sp

from alggeom.components import associated_primes, primary_decomposition, radical
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme

x, y = sp.symbols("x y")
ring = PolynomialRing((x, y))
ideal = ring.ideal(x * x, x * y)
components = primary_decomposition(ideal)
assert reduce(lambda a, b: a.intersection(b), components) == ideal
assert len(components) == 2
assert radical(ideal) == ring.ideal(x)
primes = associated_primes(ideal)
assert len(primes) == 2
assert any(p == ring.ideal(x, y) for p in primes)
assert tuple(s.ideal for s in AffineScheme(ideal).irreducible_components()) == (ring.ideal(x),)
print("Primary components:", [q.generators for q in components])
print("Associated primes:", [p.generators for p in primes])
