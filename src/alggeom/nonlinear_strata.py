"""Exact parameter strata for nonlinear polynomial hypersurface fibers.

A univariate fiber is stratified by its specialized polynomial degree and
Sylvester-matrix rank, which determines the number of distinct geometric roots.
Multivariate hypersurfaces are stratified by coefficient vanishing: a nonzero
nonconstant hypersurface has dimension n-1 over an algebraic closure.
"""

from dataclasses import dataclass
from math import comb

import sympy as sp

from .constructible import ConstructibleSet, LocallyClosedSet
from .families import FiberInvariants
from .ideals import PolynomialRing
from .linear_strata import _minors


@dataclass(frozen=True)
class NonlinearFiberPiece:
    """A constructible parameter locus with fixed geometric fiber invariants."""

    locus: ConstructibleSet
    invariants: FiberInvariants

    @property
    def dimension(self):
        """Geometric dimension of a fiber on this parameter locus."""
        return self.invariants.dimension

    @property
    def geometric_points(self):
        """Number of distinct geometric points, or infinity for positive dimension."""
        return self.invariants.geometric_points


@dataclass(frozen=True)
class NonlinearFiberStratification:
    """A complete coefficient-rank partition for one polynomial equation over QQ."""

    equation: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    pieces: tuple[NonlinearFiberPiece, ...]

    def specialize(self, values):
        """Evaluate exact geometric dimension and distinct-point cardinality."""
        if isinstance(values, (int, sp.Expr)) and len(self.parameters) == 1:
            values = (values,)
        values = tuple(values)
        if len(values) != len(self.parameters):
            raise ValueError("one value is required for every parameter")
        substitution = dict(zip(self.parameters, values, strict=True))
        polynomial = sp.Poly(
            self.equation.subs(substitution, simultaneous=True), *self.variables, extension=True
        )
        if polynomial.is_zero:
            return FiberInvariants(len(self.variables), sp.oo)
        if polynomial.total_degree() == 0:
            return FiberInvariants(-1, 0)
        if len(self.variables) > 1:
            return FiberInvariants(len(self.variables) - 1, sp.oo)
        return FiberInvariants(0, int(polynomial.sqf_part().degree()))


def _add_piece(pieces, locus, invariants):
    if locus.pieces and any(not piece.closure_ideal.is_unit for piece in locus.pieces):
        pieces.append(NonlinearFiberPiece(locus, invariants))


def _sylvester_matrix(polynomial, coordinate, degree):
    """Build the coefficient matrix of p*A + p'*B at the fixed degree."""
    derivative = sp.diff(polynomial, coordinate)
    rows = []
    width = 2 * degree - 1
    for shift in range(degree - 1):
        p = sp.Poly(polynomial * coordinate**shift, coordinate)
        rows.append([p.nth(i) for i in range(width)])
    for shift in range(degree):
        p = sp.Poly(derivative * coordinate**shift, coordinate)
        rows.append([p.nth(i) for i in range(width)])
    return sp.Matrix(rows)


def hypersurface_fiber_stratification(equation, variables, parameters, *, max_minors=4096):
    """Stratify arbitrary polynomial hypersurface fibers over QQ exactly.

    Fiber equations may be nonlinear in any number of variables and contain
    any number of polynomial parameters. In one fiber variable, Sylvester
    ranks yield exact counts of distinct geometric points, including multiple
    roots and degree drops. With at least two fiber coordinates, the
    hypersurface dimension is exact; positive-dimensional fibers have
    infinitely many geometric points. The determinant budget is explicit.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    if not variables or len(set((*variables, *parameters))) != len(variables) + len(parameters):
        raise ValueError("fiber variables and parameters must be distinct")
    expression = sp.sympify(equation)
    if expression.has(sp.Float) or expression.free_symbols - set((*variables, *parameters)):
        raise ValueError("an exact polynomial in the declared variables is required")
    polynomial = sp.Poly(expression, *variables, *parameters, domain=sp.QQ)
    if polynomial.total_degree() < 0:
        ring = PolynomialRing(parameters)
        return NonlinearFiberStratification(
            expression,
            variables,
            parameters,
            (
                NonlinearFiberPiece(
                    ConstructibleSet((LocallyClosedSet(ring.ideal()),)),
                    FiberInvariants(len(variables), sp.oo),
                ),
            ),
        )
    ring = PolynomialRing(parameters)
    fibers = sp.Poly(expression, *variables)
    max_degree = max(sum(monomial) for monomial, _ in fibers.terms())
    by_degree = {}
    for monomial, coefficient in fibers.terms():
        by_degree.setdefault(sum(monomial), []).append(sp.expand(coefficient))
    pieces = []
    for degree in range(max_degree, -1, -1):
        higher = [
            coeff for power, coeffs in by_degree.items() if power > degree for coeff in coeffs
        ]
        equalities = ring.ideal(*higher)
        if equalities.is_unit:
            continue
        leading = by_degree.get(degree, [])
        if not leading:
            continue
        degree_loci = tuple(LocallyClosedSet(equalities, (coefficient,)) for coefficient in leading)
        if degree == 0:
            _add_piece(pieces, ConstructibleSet(degree_loci), FiberInvariants(-1, 0))
        elif len(variables) > 1:
            _add_piece(
                pieces, ConstructibleSet(degree_loci), FiberInvariants(len(variables) - 1, sp.oo)
            )
        else:
            coordinate = variables[0]
            size = 2 * degree - 1
            matrix = _sylvester_matrix(expression, coordinate, degree)
            used = 0
            for rank in range(degree - 1, size + 1):
                used += comb(size, rank) ** 2 + (comb(size, rank + 1) ** 2 if rank < size else 0)
                if used > max_minors:
                    raise ValueError("Sylvester-rank stratification exceeds max_minors")
                vanishing = _minors(matrix, rank + 1)
                minor_loci = _minors(matrix, rank)
                if rank == 0:
                    minor_loci = (sp.S.One,)
                rank_ideal = equalities + ring.ideal(*vanishing)
                conditions = tuple(
                    LocallyClosedSet(rank_ideal, (leading_coeff, minor))
                    for leading_coeff in leading
                    for minor in minor_loci
                )
                _add_piece(
                    pieces,
                    ConstructibleSet(conditions),
                    FiberInvariants(0, rank - degree + 1),
                )
    # The identically-zero fiber polynomial is the common zero locus of all
    # polynomial coefficients in the fiber coordinates.
    all_coefficients = [coefficient for group in by_degree.values() for coefficient in group]
    zero = ring.ideal(*all_coefficients)
    _add_piece(
        pieces,
        ConstructibleSet((LocallyClosedSet(zero),)),
        FiberInvariants(len(variables), sp.oo),
    )
    return NonlinearFiberStratification(expression, variables, parameters, tuple(pieces))


__all__ = [
    "NonlinearFiberPiece",
    "NonlinearFiberStratification",
    "hypersurface_fiber_stratification",
]


@dataclass(frozen=True)
class SeparableFiberStratification:
    """Exact product stratification of uncoupled nonlinear fiber equations."""

    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    factors: tuple[NonlinearFiberStratification, ...]
    pieces: tuple[NonlinearFiberPiece, ...]

    def specialize(self, values):
        """Combine exact invariants of the independent univariate fibers."""
        invariants = [factor.specialize(values) for factor in self.factors]
        if any(result.dimension == -1 for result in invariants):
            return FiberInvariants(-1, 0)
        dimension = sum(result.dimension for result in invariants)
        count = sp.oo if dimension else sp.prod(result.geometric_points for result in invariants)
        return FiberInvariants(dimension, count)


def separable_fiber_stratification(equations, variables, parameters, *, max_pieces=4096):
    """Stratify products of independent nonlinear univariate equations.

    Each fiber coordinate may appear in at most one defining equation, and
    each equation may depend on at most one fiber coordinate. Coefficients can
    be arbitrary polynomials in all declared parameters. This includes
    nonreduced and singular zero-dimensional fibers. Coupled systems require
    a comprehensive multivariate Gröbner algorithm and are not accepted.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    equations = tuple(sp.sympify(e) for e in equations)
    if not variables or len(set((*variables, *parameters))) != len(variables) + len(parameters):
        raise ValueError("fiber and parameter symbols must be distinct")
    if any(e.has(sp.Float) or e.free_symbols - set((*variables, *parameters)) for e in equations):
        raise ValueError("an exact polynomial in the declared variables is required")
    assigned = {}
    for equation in equations:
        active = tuple(v for v in variables if equation.has(v))
        if len(active) != 1 or active[0] in assigned:
            raise NotImplementedError("independent univariate equations are required")
        assigned[active[0]] = equation
    factors = tuple(
        hypersurface_fiber_stratification(
            assigned.get(variable, sp.S.Zero), (variable,), parameters
        )
        for variable in variables
    )
    ring = PolynomialRing(parameters)
    states = [(ConstructibleSet((LocallyClosedSet(ring.ideal()),)), FiberInvariants(0, 1))]
    for factor in factors:
        next_states = []
        for locus, accumulated in states:
            for piece in factor.pieces:
                intersection = locus & piece.locus
                if not any(not p.closure_ideal.is_unit for p in intersection.pieces):
                    continue
                if accumulated.dimension < 0 or piece.dimension < 0:
                    invariants = FiberInvariants(-1, 0)
                else:
                    dimension = accumulated.dimension + piece.dimension
                    count = (
                        sp.oo
                        if dimension
                        else accumulated.geometric_points * piece.geometric_points
                    )
                    invariants = FiberInvariants(dimension, count)
                next_states.append((intersection, invariants))
                if len(next_states) > max_pieces:
                    raise ValueError("separable fiber stratification exceeds max_pieces")
        states = next_states
    return SeparableFiberStratification(
        equations,
        variables,
        parameters,
        factors,
        tuple(NonlinearFiberPiece(locus, inv) for locus, inv in states),
    )


__all__ += ["SeparableFiberStratification", "separable_fiber_stratification"]
