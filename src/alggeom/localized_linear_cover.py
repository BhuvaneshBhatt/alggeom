"""Complete, certified parameter covers for linear univariate fiber systems.

Each fiber equation has the form f_i(a) * x + g_i(a), over an arbitrary
QQ[a]/I and any explicitly localized basic open D(h).  The recursion splits
on coefficient ideals and residual determinants without requiring I prime,
reduced, equidimensional, or free from embedded components.  It terminates
in at most two coefficient branches per equation plus a disjoint partition
of the determinant exceptional ideal, subject to a declared piece budget.

This is a *pointwise* comprehensive cover on geometric parameter points.
On nonreduced parameter schemes, equality of geometric fibers does not imply
flat or scheme-theoretic constancy across nilpotent infinitesimal directions.
Higher-degree fiber polynomials are not covered by this module.
"""

from dataclasses import dataclass

import sympy as sp

from .parameter_localization import ParameterLocalization


class LinearCoverLimitError(RuntimeError):
    """The exact branch or verification budget was exceeded."""


def _locus(localization):
    from .constructible import ConstructibleSet, LocallyClosedSet

    return ConstructibleSet((LocallyClosedSet(localization.ideal, localization.inverted),))


def _empty(locus):
    # Saturation, unlike an approximate point search, tests emptiness exactly.
    return all(piece.ideal.saturate(sp.prod(piece.nonzero)).is_unit for piece in locus.pieces)


@dataclass(frozen=True)
class LinearCoverPiece:
    """A disjoint geometric parameter locus and one certified fiber basis."""

    localization: ParameterLocalization
    kind: str
    pivot: int | None = None
    second: int | None = None

    @property
    def condition(self):
        """Exact polynomial equalities and nonvanishing guards for this piece."""
        return self.localization.condition


@dataclass(frozen=True)
class LocalizedLinearCover:
    """Recursive complete cover for finitely many degree <= 1 equations in x."""

    equations: tuple[sp.Expr, ...]
    variable: sp.Symbol
    base: ParameterLocalization
    pieces: tuple[LinearCoverPiece, ...]

    def basis(self, piece):
        """Provide (1), (), or a monic linear polynomial over the localization."""
        if piece not in self.pieces:
            raise ValueError("unknown cover piece")
        if piece.kind == "unit":
            return (sp.S.One,)
        if piece.kind == "zero":
            return ()
        if piece.kind != "pivot" or piece.pivot is None:
            raise ArithmeticError("invalid basis piece")
        f, g = _coefficients(self.equations, self.variable, self.base)[piece.pivot]
        return (sp.cancel(self.variable + g / f),)

    def verify(self, *, max_pairs=400):
        """Independently replay every fiber identity and cover partition.

        No candidate Groebner backend is invoked. Nonvanishing witnesses are
        explicit inverted parameters; every vanishing identity is reduced
        modulo the exact saturation I:S^infinity. Coverage and disjointness
        use symbolic constructible-set differences and ideal saturation.
        """
        try:
            if not isinstance(max_pairs, int) or max_pairs < 1:
                return False
            coeffs = _coefficients(self.equations, self.variable, self.base)
            base = _locus(self.base)
            loci = []
            for piece in self.pieces:
                loc = piece.localization
                if loc.ideal.ring != self.base.ideal.ring:
                    return False
                current = _locus(loc)
                if not _empty(current - base):
                    return False
                if piece.kind == "zero":
                    if any(loc.reduce(f) != 0 or loc.reduce(g) != 0 for f, g in coeffs):
                        return False
                    if piece.pivot is not None or piece.second is not None:
                        return False
                elif piece.kind == "pivot":
                    if piece.pivot is None or not 0 <= piece.pivot < len(coeffs):
                        return False
                    if piece.second is not None:
                        return False
                    f, g = coeffs[piece.pivot]
                    if not any(loc.reduce(f - h) == 0 for h in loc.inverted):
                        return False
                    if any(
                        loc.reduce(f * other_g - other_f * g) != 0 for other_f, other_g in coeffs
                    ):
                        return False
                elif piece.kind == "unit":
                    if piece.pivot is None or not 0 <= piece.pivot < len(coeffs):
                        return False
                    f, g = coeffs[piece.pivot]
                    if piece.second is None:
                        # A constant source equation g_i is invertible.
                        if loc.reduce(f) != 0 or not any(
                            loc.reduce(g - h) == 0 for h in loc.inverted
                        ):
                            return False
                    else:
                        if not 0 <= piece.second < len(coeffs) or piece.second == piece.pivot:
                            return False
                        other_f, other_g = coeffs[piece.second]
                        determinant = f * other_g - other_f * g
                        if not any(loc.reduce(determinant - h) == 0 for h in loc.inverted):
                            return False
                else:
                    return False
                for previous in loci:
                    if not _empty(current & previous):
                        return False
                loci.append(current)
                if len(loci) > max_pairs:
                    raise LinearCoverLimitError("pairwise proof budget exceeded")
            remaining = base
            for locus in loci:
                remaining = remaining - locus
            return _empty(remaining)
        except (
            NotImplementedError,
            ValueError,
            TypeError,
            ArithmeticError,
            ZeroDivisionError,
            sp.PolynomialError,
        ):
            return False

    def specialize(self, values):
        """Select the exact fiber basis for an algebraic parameter point."""
        values = tuple(sp.sympify(value) for value in values)
        if len(values) != len(self.base.ideal.ring.variables) or any(
            value.has(sp.Float) or value.is_algebraic is not True for value in values
        ):
            raise ValueError("one exact algebraic parameter value per coordinate is required")
        extension = (
            sp.QQ.algebraic_field(*values)
            if any(value.is_Rational is not True for value in values)
            else sp.QQ
        )
        matched = [
            piece
            for piece in self.pieces
            if _locus(piece.localization).contains(values, field=extension)
        ]
        if len(matched) != 1:
            raise ValueError("not exactly one cover piece contains this parameter point")
        subs = dict(zip(self.base.ideal.ring.variables, values, strict=True))
        return tuple(sp.cancel(eq.subs(subs, simultaneous=True)) for eq in self.basis(matched[0]))


def _coefficients(equations, variable, base):
    params = base.ideal.ring.variables
    answer = []
    for eq in equations:
        poly = sp.Poly(eq, variable, *params, domain=sp.QQ)
        if poly.degree(variable) > 1:
            raise NotImplementedError("only fiber degree <= 1 is supported")
        f = sp.expand(sp.diff(eq, variable))
        g = sp.expand(eq.subs(variable, 0))
        if variable in f.free_symbols or variable in g.free_symbols:
            raise ValueError("the fiber equations must be linear")
        answer.append((f, g))
    return tuple(answer)


def comprehensive_localized_linear_cover(equations, variable, base, *, max_pieces=128):
    """Complete recursion for linear fiber equations over arbitrary QQ[a]/I.

    Branches first on the coefficient f_i, then on its constant g_i. If f_i
    is invertible, the remaining equations reduce to determinant conditions
    f_i*g_j-f_j*g_i. Splitting their *entire arbitrary ideal* yields disjoint
    unit branches and one pivot branch. Every leaf is independently verified.
    """
    if not isinstance(base, ParameterLocalization):
        raise TypeError("base must be a parameter localization")
    if not isinstance(variable, sp.Symbol) or variable in base.ideal.ring.variables:
        raise ValueError("fiber variable must be distinct from parameters")
    if not isinstance(max_pieces, int) or max_pieces < 1:
        raise ValueError("max_pieces must be positive")
    equations = tuple(sp.expand(sp.sympify(eq)) for eq in equations)
    coefficients = _coefficients(equations, variable, base)
    pieces = []

    def emit(localization, kind, pivot=None, second=None):
        if len(pieces) >= max_pieces:
            raise LinearCoverLimitError("parameter piece budget exhausted")
        pieces.append(LinearCoverPiece(localization, kind, pivot, second))

    def branch(localization, index):
        if index == len(coefficients):
            emit(localization, "zero")
            return
        f, g = coefficients[index]
        open_f, closed_f = localization.split_at(f) if f != 0 else (None, localization)
        if open_f is not None:
            residuals = tuple(
                (j, sp.expand(f * gj - fj * g))
                for j, (fj, gj) in enumerate(coefficients)
                if j != index
            )
            # A disjoint ordered principal-open decomposition of D(J),
            # J=(f_i*g_j-f_j*g_i). We retain the exact source determinant
            # for every unit witness; canonical ideal generators may differ
            # from the original determinant list, so never infer a witness
            # from the ideal's internally normalized generators.
            pivot_branch = open_f
            for j, determinant in residuals:
                if pivot_branch is None:
                    break
                unit_branch, pivot_branch = pivot_branch.split_at(determinant)
                if unit_branch is not None:
                    emit(unit_branch, "unit", index, j)
            if pivot_branch is not None:
                emit(pivot_branch, "pivot", index)
        if closed_f is not None:
            open_g, closed_g = closed_f.split_at(g) if g != 0 else (None, closed_f)
            if open_g is not None:
                emit(open_g, "unit", index)
            if closed_g is not None:
                branch(closed_g, index + 1)

    branch(base, 0)
    cover = LocalizedLinearCover(equations, variable, base, tuple(pieces))
    if not cover.verify(max_pairs=max(max_pieces, 1)):
        raise ArithmeticError("localized linear cover failed independent verification")
    return cover


__all__ = [
    "LinearCoverLimitError",
    "LinearCoverPiece",
    "LocalizedLinearCover",
    "comprehensive_localized_linear_cover",
]
