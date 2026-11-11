"""Reconstruct implicit compositions and check resultant multiplicities."""

import sympy as sp
from algroots.quotient import QuotientAlgebra

from alggeom.decomposition.algebraic import reconstruct_algebraic

x, z, a, t = sp.symbols("x z a t")
relation = x**4 - z * x**2 * (x + 1) + z**3 * (x + 1) ** 2
ambient = QuotientAlgebra.from_polynomials([a**2 - z * a + z**3], (a,), domain=sp.QQ.frac_field(z))
result = reconstruct_algebraic(
    relation, (x, z), candidate=x**2 - a * (x + 1), ambient=ambient, parameter=t
)
assert result.verify() and not result.coefficient_field.contains(z)
decomposition = result.decomposition
assert decomposition.outer == sp.Poly(x**2 + t * (x + 1), x, t, domain=sp.QQ)
assert decomposition.inner == sp.Poly(t**2 + z * t + z**3, t, z, domain=sp.QQ)
assert (
    sp.expand(
        sp.resultant(decomposition.outer.as_expr(), decomposition.inner.as_expr(), t) - relation
    )
    == 0
)

# Each root of the cubic belongs to two of its three conjugate quadratic factors.
cubic = x**3 - 3 * x - z
cubic_ambient = QuotientAlgebra.from_polynomials(
    [a**3 - 3 * a - z], (a,), domain=sp.QQ.frac_field(z)
)
power_result = reconstruct_algebraic(
    cubic, (x, z), candidate=x**2 + a * x + a**2 - 3, ambient=cubic_ambient, parameter=t
)
assert power_result.verify() and power_result.decomposition.exponent == 2
assert sp.expand(sp.resultant(x**2 + t * x + t**2 - 3, t**3 - 3 * t - z, t) - cubic**2) == 0
print("Coefficient field excludes the input:", not result.coefficient_field.contains(z))
print("Cubic resultant power:", power_result.decomposition.exponent)
