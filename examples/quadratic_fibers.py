"""Separate generic finite algebra from degree drops in a family of quadratic fibers."""

import sympy as sp
from algroots.quotient import QuotientAlgebra
from semialg.decision import equivalent
from semialg.qe import quantifier_eliminate

from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme

x, a, b, c = sp.symbols("x a b c", real=True)
f = a * x * x + b * x + c
generic = QuotientAlgebra.from_polynomials((f,), (x,), domain=sp.QQ.frac_field(a, b, c))
assert generic.dimension == 2
S = PolynomialRing((x,))
assert (
    AffineScheme(S.ideal(f.subs({a: 0, b: 2, c: -1}))).coordinate_ring.vector_space_dimension == 1
)
assert AffineScheme(S.ideal(f.subs({a: 0, b: 0, c: 0}))).coordinate_ring.krull_dimension == 1
assert AffineScheme(S.ideal(f.subs({a: 0, b: 0, c: 1}))).is_empty
u, d = sp.symbols("u d", real=True)
radicand = quantifier_eliminate(sp.Eq(u * u, d), [("exists", u)], variables=(d, u))
assert equivalent(radicand, d >= 0, (d,))
assert sp.expand(4 * a * f - ((2 * a * x + b) ** 2 - (b * b - 4 * a * c))) == 0
condition = sp.Or(
    sp.And(sp.Ne(a, 0), radicand.subs(d, b * b - 4 * a * c)),
    sp.And(sp.Eq(a, 0), sp.Ne(b, 0)),
    sp.And(sp.Eq(a, 0), sp.Eq(b, 0), sp.Eq(c, 0)),
)
assert bool(condition.subs({a: 1, b: 0, c: -1}))
assert not bool(condition.subs({a: 1, b: 0, c: 1}))
assert bool(condition.subs({a: 0, b: 2, c: -1}))
print(
    "Generic fibers have length two; linear, empty and entire-line fibers need separate parameter strata."
)
