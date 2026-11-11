"""Exact geometric fiber stratification via constructive algebraic elimination.

Every input is a finite polynomial equation/inequation system over QQ, with
arbitrary coupling among the fiber coordinates. The algorithm computes the
fiber-dimension filtration by detecting dense coordinate projections. It can
also partition the zero-dimensional locus by its number of geometric points,
using an explicitly bounded sequence of existential distinct-point tests.

These are complete constructible conditions whenever the branch and cardinality
budgets allow completion. This is a quantifier-elimination algorithm, not a
comprehensive Groebner basis: it does not return parametric Groebner bases.
"""

from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from .chevalley import project_constructible
from .constructible import ConstructibleSet, LocallyClosedSet
from .families import FiberInvariants
from .ideals import PolynomialRing


class StratificationLimitError(RuntimeError):
    """The requested complete stratification exceeds a declared degree/branch limit."""


def _union(pieces):
    return ConstructibleSet(tuple(piece for group in pieces for piece in group.pieces))


def _coefficient_conditions(expression, variables):
    """Return coefficients whose simultaneous vanishing makes a polynomial zero."""
    if not variables:
        return (sp.expand(expression),)
    return tuple(sp.Poly(expression, *variables).coeffs())


def _dense_parameter_locus(image, projections, parameters):
    """Parameters where a constructible subset of A^d is Zariski dense.

    For V(I) intersect D(g_1)...D(g_r), the specialized fiber is dense in
    irreducible affine space exactly when every equation specializes to the
    zero polynomial and no guard specializes to the zero polynomial.
    """
    ring = PolynomialRing(parameters)
    results = []
    for piece in image.pieces:
        coefficients = [
            coefficient
            for eq in piece.ideal.generators
            for coefficient in _coefficient_conditions(eq, projections)
        ]
        ideal = ring.ideal(*coefficients)
        options = [()]  # Distribute the nonzero-coefficient choices for each guard.
        for guard in piece.nonzero:
            coeffs = _coefficient_conditions(guard, projections)
            options = [(*choice, coefficient) for choice in options for coefficient in coeffs]
        for guards in options:
            if not ideal.saturate(sp.prod(guards)).is_unit:
                results.append(LocallyClosedSet(ideal, guards))
    return ConstructibleSet(tuple(results))


def _projection_source(equations, variables, parameters, nonzero):
    ring = PolynomialRing((*variables, *parameters))
    return ConstructibleSet((LocallyClosedSet(ring.ideal(*equations), nonzero),))


def _eliminate_linear_coordinates(equations, variables, guards):
    """Eliminate a coordinate with a nonzero constant linear coefficient.

    Each step is an isomorphism between the original solution set and its
    reduced graph, so it preserves cardinality and distinctness.
    """
    equations = list(equations)
    guards = list(guards)
    coordinates = list(variables)
    reconstructed = {symbol: symbol for symbol in variables}
    while True:
        chosen = None
        for index, equation in enumerate(equations):
            for variable in coordinates:
                derivative = sp.diff(equation, variable)
                if derivative.is_Rational and derivative != 0:
                    chosen = (index, variable, sp.expand(-equation.subs(variable, 0) / derivative))
                    break
            if chosen:
                break
        if chosen is None:
            break
        index, variable, replacement = chosen
        equations = [
            sp.expand(eq.subs(variable, replacement))
            for position, eq in enumerate(equations)
            if position != index
        ]
        guards = [sp.expand(guard.subs(variable, replacement)) for guard in guards]
        reconstructed = {
            symbol: sp.expand(expr.subs(variable, replacement))
            for symbol, expr in reconstructed.items()
        }
        coordinates.remove(variable)
    return tuple(equations), tuple(coordinates), tuple(guards), reconstructed


def _exists_at_least(equations, variables, parameters, nonzero, count, max_branches):
    """Decide existence of k pairwise different points using a separator.

    Over an infinite field there exists a linear functional separating any
    finite set of pairwise distinct points. This replaces a disjunction of
    coordinate inequalities for each pair with one polynomial inequation.
    """
    if count == 1:
        source = _projection_source(equations, variables, parameters, nonzero)
        return project_constructible(source, variables, max_branches=max_branches)
    reduced_equations, reduced_variables, reduced_guards, reconstruction = (
        _eliminate_linear_coordinates(equations, variables, nonzero)
    )
    copied = tuple(
        tuple(sp.Dummy(f"point_{index}_{j}") for j in range(len(reduced_variables)))
        for index in range(count)
    )
    equations_copied = []
    guards_copied = []
    for point in copied:
        replacements = dict(zip(reduced_variables, point, strict=True))
        equations_copied.extend(expr.xreplace(replacements) for expr in reduced_equations)
        guards_copied.extend(expr.xreplace(replacements) for expr in reduced_guards)
    hidden = tuple(symbol for point in copied for symbol in point)
    ring = PolynomialRing((*hidden, *parameters))
    result = []
    # A nonzero polynomial of degree at most n-1 has at most n-1 integer
    # roots. For each distinct pair, only that many linear forms
    # L_t(x)=sum_j t**j*x_j fail to separate it. Pigeonhole gives the
    # following finite collection a separating form for every k-point set.
    attempts = (len(variables) - 1) * count * (count - 1) // 2 + 1
    for shift in range(attempts):
        values = [
            sp.Add(
                *(
                    sp.Integer(shift) ** index
                    * reconstruction[variable].xreplace(
                        dict(zip(reduced_variables, point, strict=True))
                    )
                    for index, variable in enumerate(variables)
                )
            )
            for point in copied
        ]
        separator = sp.prod(
            values[i] - values[j] for i in range(count) for j in range(i + 1, count)
        )
        source = ConstructibleSet(
            (LocallyClosedSet(ring.ideal(*equations_copied), (*guards_copied, separator)),)
        )
        result.append(project_constructible(source, hidden, max_branches=max_branches))
    return _union(result)


@dataclass(frozen=True)
class GeometricFiberPiece:
    """A constructible parameter locus with fixed dimension and point count."""

    locus: ConstructibleSet
    invariants: FiberInvariants


@dataclass(frozen=True)
class ComprehensiveFiberStratification:
    """Exact constructible dimension and cardinality partition of affine fibers.

    All geometric quantities refer to the algebraic closure of QQ. The pieces
    may be represented by overlapping locally closed terms internally, but
    distinct invariant values occur on disjoint parameter loci.
    """

    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    nonzero: tuple[sp.Expr, ...]
    pieces: tuple[GeometricFiberPiece, ...]

    def specialize(self, values):
        """Return the certified invariants of an exact algebraic parameter tuple."""
        if isinstance(values, (int, sp.Expr)) and len(self.parameters) == 1:
            values = (values,)
        values = tuple(map(sp.sympify, values))
        if len(values) != len(self.parameters) or any(v.is_algebraic is not True for v in values):
            raise ValueError("an exact algebraic value is required for each parameter")
        field = (
            sp.QQ
            if all(value.is_Rational is True for value in values)
            else sp.QQ.algebraic_field(*values)
        )
        match = [
            piece.invariants for piece in self.pieces if piece.locus.contains(values, field=field)
        ]
        if not match or any(item != match[0] for item in match):
            raise ArithmeticError("parameter specialization did not match a unique invariant")
        return match[0]


def comprehensive_fiber_stratification(
    equations,
    variables,
    parameters,
    *,
    nonzero=(),
    max_branches=4096,
    count_points=False,
    max_points=4,
):
    """Partition arbitrary coupled polynomial families by exact fiber invariants.

    Constructs the dimension filtration using Chevalley elimination for every
    coordinate projection. Over an algebraically closed field, a fiber has
    dimension at least d iff one projection onto d of its coordinates is dense.
    For zero-dimensional fibers, at least-k distinct-point loci are computed
    by existential elimination and refined into exact cardinality strata.

    `max_points` bounds the Bézout upper bound, not a sampling or guessing
    cutoff. When a requested complete answer exceeds this bound or Chevalley
    exceeds `max_branches`, the function raises without returning partial data.
    Set `count_points=False` to compute dimensions alone for expensive inputs.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    if len(set((*variables, *parameters))) != len(variables) + len(parameters) or not variables:
        raise ValueError("nonempty, disjoint fiber and parameter coordinates are required")
    expressions = tuple(map(sp.sympify, equations))
    guards = tuple(map(sp.sympify, nonzero))
    if any(
        expr.has(sp.Float) or expr.free_symbols - set((*variables, *parameters))
        for expr in (*expressions, *guards)
    ):
        raise ValueError("exact polynomials in the declared coordinates are required")
    try:
        for expr in (*expressions, *guards):
            sp.Poly(expr, *variables, *parameters, domain=sp.QQ)
    except (sp.PolynomialError, sp.CoercionFailed) as exc:
        raise ValueError("polynomial equations and inequations over QQ are required") from exc
    if not isinstance(max_branches, int) or max_branches < 1:
        raise ValueError("max_branches must be positive")
    if not isinstance(max_points, int) or max_points < 1:
        raise ValueError("max_points must be positive")

    ambient = PolynomialRing(parameters)
    full = ConstructibleSet((LocallyClosedSet(ambient.ideal()),))
    source = _projection_source(expressions, variables, parameters, guards)
    filtration = {}
    for dimension in range(len(variables) + 1):
        projected = []
        for coordinate_set in combinations(variables, dimension):
            forgotten = tuple(v for v in variables if v not in coordinate_set)
            image = project_constructible(source, forgotten, max_branches=max_branches)
            projected.append(_dense_parameter_locus(image, coordinate_set, parameters))
        filtration[dimension] = _union(projected)

    pieces = [GeometricFiberPiece(full - filtration[0], FiberInvariants(-1, 0))]
    for dimension in range(1, len(variables) + 1):
        exactly = filtration[dimension] - (
            filtration[dimension + 1] if dimension < len(variables) else ConstructibleSet(())
        )
        pieces.append(GeometricFiberPiece(exactly, FiberInvariants(dimension, sp.oo)))

    zero_dimensional = filtration[0] - filtration[1]
    if not count_points:
        pieces.append(GeometricFiberPiece(zero_dimensional, FiberInvariants(0, None)))
    else:
        degrees = [sp.Poly(eq, *variables).total_degree() for eq in expressions]
        upper_bound = max(
            1, max((max(1, degree) for degree in degrees), default=1) ** len(variables)
        )
        if upper_bound > max_points:
            raise StratificationLimitError(
                f"zero-dimensional Bézout bound {upper_bound} exceeds max_points={max_points}"
            )
        at_least = {1: filtration[0]}
        for count in range(2, upper_bound + 1):
            at_least[count] = _exists_at_least(
                expressions, variables, parameters, guards, count, max_branches
            )
        for count in range(1, upper_bound + 1):
            count_locus = at_least[count] - (
                at_least[count + 1] if count < upper_bound else ConstructibleSet(())
            )
            pieces.append(
                GeometricFiberPiece(zero_dimensional & count_locus, FiberInvariants(0, count))
            )
    return ComprehensiveFiberStratification(
        expressions, variables, parameters, guards, tuple(pieces)
    )


__all__ = [
    "StratificationLimitError",
    "GeometricFiberPiece",
    "ComprehensiveFiberStratification",
    "comprehensive_fiber_stratification",
]
