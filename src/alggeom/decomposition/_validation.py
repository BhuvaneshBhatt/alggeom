"""Exact coefficient domains, symbols, relations and visible pole guards."""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import AlgebraicField
from sympy.polys.polyerrors import BasePolynomialError

from .._algebra import visible_denominators as visible_denominators


class DecompositionError(ValueError):
    """Invalid decomposition input or a failed exact reconstruction."""


class DegreeLimitError(DecompositionError):
    """The input exceeds a pre-computation degree limit."""


def validate_symbols(variables):
    if any(not isinstance(v, sp.Symbol) for v in variables) or len(set(variables)) != len(
        variables
    ):
        raise DecompositionError("distinct Symbol variables are required")


def validate_integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise DecompositionError(f"{name} must be an integer >= {minimum}")


def normalize_relation(expression, variables, constants):
    expression = sp.sympify(expression)
    if expression.has(sp.Float) or expression.free_symbols - set(variables):
        raise DecompositionError("an exact relation over the constants is required")
    numerator = sp.fraction(sp.cancel(expression))[0]
    return sp.Poly(numerator, *variables, domain=constants).monic()


def parse_relation(expression, variables, constants):
    """Validate an irreducible exact polynomial depending on both coordinates."""
    expression = sp.sympify(expression)
    if expression.has(sp.Float) or expression.free_symbols - set(variables):
        raise DecompositionError("an exact polynomial relation is required")
    try:
        polynomial = sp.Poly(expression, *variables, domain=constants)
    except (ValueError, TypeError, BasePolynomialError) as exc:
        raise DecompositionError("an exact polynomial relation is required") from exc
    if any(polynomial.degree(v) <= 0 for v in variables) or not polynomial.is_irreducible:
        raise DecompositionError("the relation must be irreducible and depend on both variables")
    return polynomial.monic()


def constant_field(domain):
    if domain == sp.ZZ:
        return sp.QQ
    if domain == sp.QQ or isinstance(domain, AlgebraicField):
        return domain
    raise DecompositionError("constants must be QQ or an exact algebraic number field")


def input_variables(variables):
    values = (variables,) if isinstance(variables, sp.Symbol) else tuple(variables)
    validate_symbols(values)
    if not values or any(not value.is_commutative for value in values):
        raise DecompositionError("nonempty commutative input coordinates are required")
    return values


def validate_real_symbols(symbols, name):
    values = tuple(symbols)
    validate_symbols(values)
    if any(value.is_real is False for value in values):
        raise DecompositionError(f"{name} must be real coordinates")
    return values


class IncompleteDecompositionError(DecompositionError):
    """A requested complete mathematical result exceeds the search scope.

    The partial search or chain certificate is available as ``result``.
    """

    def __init__(self, result):
        self.result = result
        super().__init__(
            "decomposition search is incomplete; inspect error.result for its evidence"
        )
