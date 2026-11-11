"""Constructive existential elimination over algebraically closed characteristic-zero fields.

The implementation uses comprehensive univariate polynomial gcds. Parameter
coefficient vanishing is split into locally closed loci, and pseudoremainders
are computed only on branches where the divisor's leading coefficient is
nonzero. Repeating the one-variable projection eliminates arbitrary finite
sequences of existential polynomial variables. All computations are exact over
QQ and interpreted geometrically over its algebraic closure.
"""

from dataclasses import dataclass

import sympy as sp

from .constructible import ConstructibleSet, LocallyClosedSet
from .ideals import PolynomialRing


class ProjectionLimitError(RuntimeError):
    """The constructive projection exceeded its explicit branch budget."""


@dataclass
class _BranchBudget:
    maximum: int
    used: int = 0

    def visit(self):
        self.used += 1
        if self.used > self.maximum:
            raise ProjectionLimitError("constructive projection exceeded max_branches")


def _nonempty(locus):
    # Saturation decides whether V(I) meets the principal open D(product).
    return not locus.ideal.saturate(sp.prod(locus.nonzero)).is_unit


def _reduce_polynomial(expression, coordinate, ideal):
    """Reduce every coefficient modulo parameter equalities before splitting."""
    poly = sp.Poly(expression, coordinate)
    return sp.expand(
        sum(
            ideal.normal_form(coefficient) * coordinate ** exponent[0]
            for exponent, coefficient in poly.terms()
        )
    )


def _degree_cases(expression, coordinate, locus, budget):
    """Split where the specialized degree drops, retaining an exact leading term."""
    budget.visit()
    expression = _reduce_polynomial(expression, coordinate, locus.ideal)
    if expression == 0:
        return [(locus, sp.S.Zero)]
    polynomial = sp.Poly(expression, coordinate)
    leading = polynomial.LC()
    # The open branch has a stable degree; on the closed branch the leading
    # term disappears and the degree-splitting process makes progress.
    open_locus = LocallyClosedSet(locus.ideal, (*locus.nonzero, leading))
    closed_locus = LocallyClosedSet(locus.ideal + locus.ideal.ring.ideal(leading), locus.nonzero)
    result = [(open_locus, expression)] if _nonempty(open_locus) else []
    if _nonempty(closed_locus):
        result.extend(_degree_cases(expression, coordinate, closed_locus, budget))
    return result


def _pseudoremainder(first, second, coordinate, parameters):
    """Pseudodivide over a polynomial coefficient ring without denominators."""
    domain = sp.QQ.poly_ring(*parameters) if parameters else sp.QQ
    lhs = sp.Poly(first, coordinate, domain=domain)
    rhs = sp.Poly(second, coordinate, domain=domain)
    return sp.expand(lhs.prem(rhs).as_expr())


def _gcd_pair(first, second, coordinate, locus, budget):
    """Enumerate exact specialized gcds up to units on locally closed loci."""
    for current_locus, first_poly in _degree_cases(first, coordinate, locus, budget):
        for next_locus, second_poly in _degree_cases(second, coordinate, current_locus, budget):
            if first_poly == 0:
                yield next_locus, second_poly
            elif second_poly == 0:
                yield next_locus, first_poly
            else:
                left = sp.Poly(first_poly, coordinate)
                right = sp.Poly(second_poly, coordinate)
                if left.degree() < right.degree():
                    first_poly, second_poly = second_poly, first_poly
                remainder = _pseudoremainder(
                    first_poly, second_poly, coordinate, next_locus.ideal.ring.variables
                )
                yield from _gcd_pair(second_poly, remainder, coordinate, next_locus, budget)


def _gcd_cases(equations, coordinate, locus, budget):
    candidates = [(locus, sp.S.Zero)]
    for equation in equations:
        candidates = [
            outcome
            for current_locus, gcd in candidates
            for outcome in _gcd_pair(gcd, equation, coordinate, current_locus, budget)
        ]
    return candidates


def project_constructible(constructible, variables, *, max_branches=4096):
    """Project a constructible set along specified coordinates exactly.

    Input pieces are finite intersections of polynomial zero and nonzero
    conditions over QQ. The result is a finite union of locally closed sets,
    equal to the existential image over an algebraic closure of QQ. Every
    parameter-specialization branch is resolved using exact polynomial gcd
    arithmetic; hitting ``max_branches`` raises instead of returning an
    incomplete image. Complexity can be very high for nonlinear systems.
    """
    if not isinstance(constructible, ConstructibleSet):
        raise TypeError("a ConstructibleSet is required")
    if constructible.ring is None:
        raise ValueError("empty constructible sets require an ambient polynomial ring")
    if constructible.ring.domain != sp.QQ:
        raise NotImplementedError(
            "constructive projection currently requires rational coefficients"
        )
    if isinstance(variables, sp.Symbol):
        variables = (variables,)
    variables = tuple(variables)
    if len(set(variables)) != len(variables) or any(
        v not in constructible.ring.variables for v in variables
    ):
        raise ValueError("projected coordinates must be distinct ambient variables")
    if not isinstance(max_branches, int) or max_branches < 1:
        raise ValueError("max_branches must be a positive integer")
    budget = _BranchBudget(max_branches)
    result = constructible
    for coordinate in variables:
        remain = tuple(v for v in result.ring.variables if v != coordinate)
        ring = PolynomialRing(remain, sp.QQ)
        pieces = []
        for piece in result.pieces:
            locus = LocallyClosedSet(ring.ideal())
            guard = sp.prod(piece.nonzero)
            for base_locus, gcd in _gcd_cases(piece.ideal.generators, coordinate, locus, budget):
                if gcd == 0:
                    # With no equations, some coordinate avoids a polynomial
                    # guard iff that polynomial does not specialize to zero.
                    for final_locus, nonzero_poly in _degree_cases(
                        guard, coordinate, base_locus, budget
                    ):
                        if nonzero_poly != 0:
                            pieces.append(final_locus)
                    continue
                degree = sp.Poly(gcd, coordinate).degree()
                if degree == 0:
                    continue
                for final_locus, shared in _gcd_pair(
                    gcd, sp.expand(guard**degree), coordinate, base_locus, budget
                ):
                    if sp.Poly(shared, coordinate).degree() < degree:
                        pieces.append(final_locus)
        # The intermediate ring is retained even when the image is empty by
        # representing it as the empty algebraic locus V(1).
        result = ConstructibleSet(tuple(pieces) if pieces else (LocallyClosedSet(ring.ideal(1)),))
    return result


__all__ = ["ProjectionLimitError", "project_constructible"]
