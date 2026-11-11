"""Exact determinantal partitions for parameter-dependent affine-linear systems."""

from dataclasses import dataclass
from itertools import combinations
from math import comb

import sympy as sp

from .constructible import ConstructibleSet, LocallyClosedSet
from .ideals import PolynomialRing


def _minors(matrix, size):
    """Enumerate nonzero minors, exploiting the zero size as a single unit minor."""
    if size == 0:
        return (sp.S.One,)
    if size > min(matrix.shape):
        return ()
    expressions = []
    for rows in combinations(range(matrix.rows), size):
        for cols in combinations(range(matrix.cols), size):
            value = sp.factor(matrix.extract(rows, cols).det())
            if value != 0 and value not in expressions:
                expressions.append(value)
    return tuple(expressions)


@dataclass(frozen=True)
class LinearFiberPiece:
    """Constructible parameter locus with fixed affine-fiber dimension."""

    locus: ConstructibleSet
    dimension: int
    geometric_points: object
    rank: int


@dataclass(frozen=True)
class LinearFiberStratification:
    """Exact rank partition of A(a)x=b(a) over an algebraically closed field."""

    parameters: tuple[sp.Symbol, ...]
    variables: tuple[sp.Symbol, ...]
    coefficient_matrix: sp.ImmutableMatrix
    constants: sp.ImmutableMatrix
    pieces: tuple[LinearFiberPiece, ...]

    def specialize(self, values):
        """Return fiber dimension and distinct point count at an exact parameter tuple."""
        from .families import FiberInvariants

        values = tuple(values)
        if len(values) != len(self.parameters):
            raise ValueError("one value is required for every parameter")
        subs = dict(zip(self.parameters, values, strict=True))
        matrix = self.coefficient_matrix.subs(subs, simultaneous=True)
        augmented = matrix.row_join(self.constants.subs(subs, simultaneous=True))
        rank = int(matrix.rank())
        if augmented.rank() != rank:
            return FiberInvariants(-1, 0)
        dimension = len(self.variables) - rank
        return FiberInvariants(dimension, 1 if dimension == 0 else sp.oo)


def linear_fiber_stratification(equations, variables, parameters, *, max_minors=4096):
    """Partition a polynomial-parameter family of affine-linear systems exactly.

    The equations must be affine linear in the fiber coordinates, but may be
    polynomial of arbitrary degree in the parameters. Each returned piece is
    a finite union of locally closed determinantal conditions. Inconsistent
    systems have empty fibers; positive-dimensional affine fibers have
    infinitely many points over the algebraic closure.
    """
    variables = tuple(variables)
    parameters = tuple(parameters)
    if not variables or len(set((*variables, *parameters))) != len(variables) + len(parameters):
        raise ValueError(
            "fiber and parameter symbols must be distinct; fiber variables are required"
        )
    expressions = tuple(sp.sympify(e) for e in equations)
    if any(e.has(sp.Float) or e.free_symbols - set((*variables, *parameters)) for e in expressions):
        raise ValueError("exact polynomials in the declared variables are required")
    for expression in expressions:
        if sp.Poly(expression, *variables).total_degree() > 1:
            raise NotImplementedError("the family must be affine-linear in its fiber variables")
    matrix, constants = (
        sp.linear_eq_to_matrix(expressions, variables)
        if expressions
        else (sp.zeros(0, len(variables)), sp.zeros(0, 1))
    )
    augmented = matrix.row_join(constants)
    ring = PolynomialRing(parameters)
    pieces = []
    budget = 0
    for rank in range(min(matrix.shape) + 1):
        minors = _minors(matrix, rank)
        if rank == 0:
            minors = (sp.S.One,)
        budget += comb(matrix.rows, rank) * comb(matrix.cols, rank) if rank else 1
        budget += comb(augmented.rows, rank + 1) * comb(augmented.cols, rank + 1)
        if budget > max_minors:
            raise ValueError("rank stratification exceeds max_minors")
        vanishing = (*_minors(matrix, rank + 1), *_minors(augmented, rank + 1))
        equalities = ring.ideal(*vanishing)
        if equalities.is_unit:
            continue
        loci = tuple(
            LocallyClosedSet(equalities, () if rank == 0 else (minor,)) for minor in minors
        )
        if not loci:
            continue
        dimension = len(variables) - rank
        pieces.append(
            LinearFiberPiece(
                ConstructibleSet(loci), dimension, 1 if dimension == 0 else sp.oo, rank
            )
        )
    full_space = ConstructibleSet((LocallyClosedSet(ring.ideal()),))
    nonempty = ConstructibleSet(())
    for piece in pieces:
        nonempty = nonempty | piece.locus
    empty = full_space - nonempty
    if empty.pieces:
        pieces.append(LinearFiberPiece(empty, -1, 0, -1))
    return LinearFiberStratification(
        parameters,
        variables,
        sp.ImmutableMatrix(matrix),
        sp.ImmutableMatrix(constants),
        tuple(pieces),
    )


__all__ = ["LinearFiberPiece", "LinearFiberStratification", "linear_fiber_stratification"]
