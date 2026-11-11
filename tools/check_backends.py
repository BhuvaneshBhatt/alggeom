"""Check that CI exercises the requested arithmetic and optional dependency profile."""

import argparse
from importlib.util import find_spec

import sympy as sp
from algroots.finite_decomposition import primary_decomposition, reduced_algebra
from algroots.finite_fields import finite_field
from algroots.quotient import QuotientAlgebra
from sympy.external.gmpy import GROUND_TYPES


def main():
    """Validate backend selection and a finite algebra over a rational function field."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("core", "real"), required=True)
    parser.add_argument("--ground-types", choices=("flint", "python"), required=True)
    args = parser.parse_args()
    if args.ground_types != GROUND_TYPES:
        raise RuntimeError(f"expected {args.ground_types} coefficients, got {GROUND_TYPES}")
    a, z, t = sp.symbols("a z t")
    algebra = QuotientAlgebra.from_polynomials([a**2 - z], (a,), domain=sp.QQ.frac_field(z))
    if algebra.minimal_polynomial(a, t) != sp.Poly(t**2 - z, t, domain=sp.QQ.frac_field(z)):
        raise RuntimeError("the quotient backend lacks function-field minimal polynomials")
    alpha = sp.sqrt(-3)
    field = sp.QQ.algebraic_field(alpha).frac_field(z)
    imaginary = QuotientAlgebra.from_polynomials([a**2 - z - alpha], (a,), domain=field)
    if imaginary.minimal_polynomial(a, t).all_coeffs() != [1, 0, -z - alpha]:
        raise RuntimeError("the quotient backend lost its algebraic coefficient embedding")
    local = QuotientAlgebra.from_polynomials([(a - z) ** 2], (a,), domain=field)
    factors = primary_decomposition(local)
    if reduced_algebra(local).dimension != 1 or [
        (q.dimension, p.dimension) for q, p in factors
    ] != [(2, 1)]:
        raise RuntimeError("the quotient backend lacks exact finite primary factors")
    finite = QuotientAlgebra.from_polynomials((a * a, z * z), (a, z), domain=sp.GF(2))
    if finite.trace_pairing != sp.zeros(4) or reduced_algebra(finite).dimension != 1:
        raise RuntimeError("the quotient backend lacks Frobenius nilradicals")
    extension = finite_field(3, 2)
    if extension.generator**9 != extension.generator:
        raise RuntimeError("the backend failed finite-extension arithmetic")
    available = find_spec("semialg") is not None
    if available != (args.profile == "real"):
        raise RuntimeError(f"semialg availability does not match the {args.profile} profile")
    if available:
        from semialg.decision import is_satisfiable

        if is_satisfiable(z**2 < 0, (z,)):
            raise RuntimeError("the real backend failed an exact decision")
    print(f"Validated {args.profile} dependencies with {GROUND_TYPES} coefficients")


if __name__ == "__main__":
    main()
