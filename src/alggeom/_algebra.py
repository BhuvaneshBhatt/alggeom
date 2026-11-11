"""Shared exact polynomial arithmetic and coordinate validation."""

import sympy as sp
from algroots.quotient import _require_quotient_domain

exact_field = _require_quotient_domain


def coordinates(values, *, allow_empty=True):
    values = (values,) if isinstance(values, sp.Symbol) else tuple(values)
    if (
        (not values and not allow_empty)
        or len(set(values)) != len(values)
        or any(not isinstance(v, sp.Symbol) or v.is_commutative is not True for v in values)
    ):
        raise ValueError("distinct commutative Symbol coordinates are required")
    return values


def integer(value, name, minimum=0):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, sp.Integer))
        or (minimum is not None and value < minimum)
    ):
        constraint = "" if minimum is None else f" >= {minimum}"
        raise ValueError(f"{name} must be an integer{constraint}")
    return int(value)


def visible_denominators(expression):
    """Find affine denominator exclusions before rational cancellation."""
    expression = sp.sympify(expression)
    denominators = [sp.together(expression).as_numer_denom()[1]]
    for term in sp.preorder_traversal(expression):
        if term.is_Pow and term.exp.is_Integer and term.exp < 0:
            denominators.extend(sp.together(term.base).as_numer_denom())
    return tuple(dict.fromkeys(d for d in denominators if d != 1))


def rational_parts(expression, ring):
    """Reduce a rational expression over the declared coordinate and coefficient fields."""
    expression = expression.as_expr() if isinstance(expression, sp.Poly) else sp.sympify(expression)
    if expression.has(sp.Float):
        raise ValueError("exact rational functions are required")
    numerator, denominator = sp.together(expression).as_numer_denom()
    left, right = ring._poly(numerator), ring._poly(denominator)
    if not right:
        raise ZeroDivisionError("a rational function has a zero denominator")
    common = left.gcd(right)
    left, right = left.exquo(common), right.exquo(common)
    scale = right.LC
    domain = right.ring.domain
    return left.mul_ground(domain.exquo(domain.one, scale)), right.monic()


__all__ = []
