"""Locally closed algebraic sets and exact existential images.

ConstructibleSet stores a finite union of V(I) intersect D(g_1)...D(g_r).
ConstructibleImage retains an existential presentation when a finite locally
closed decomposition has not been computed. Point membership is over an
algebraic closure of the declared exact coefficient field, not merely the
rational points of that field.
"""

from dataclasses import dataclass

import sympy as sp

from .ideals import Ideal


@dataclass(frozen=True)
class LocallyClosedSet:
    """The algebraic set V(ideal) with each supplied polynomial nonzero."""

    ideal: Ideal
    nonzero: tuple[sp.Expr, ...] = ()

    def __post_init__(self):
        if not isinstance(self.ideal, Ideal):
            raise TypeError("a polynomial ideal is required")
        object.__setattr__(self, "nonzero", tuple(self.ideal.ring.convert(g) for g in self.nonzero))

    @property
    def closure_ideal(self):
        """The Zariski closure, including removal of locally empty components."""
        return self.ideal.saturate(sp.prod(self.nonzero)).radical()

    def contains(self, values, *, field=None):
        """Test an exact field-valued point, using a specified extension if needed."""
        from algroots.domains import field_element

        ring = self.ideal.ring
        domain = ring.domain if field is None else field
        coords = tuple(values)
        if len(coords) != len(ring.variables):
            raise ValueError("one point coordinate is required per variable")
        substitutions = dict(zip(ring.variables, coords, strict=True))
        return all(
            field_element(f.subs(substitutions, simultaneous=True), domain) == domain.zero
            for f in self.ideal.generators
        ) and all(
            field_element(g.subs(substitutions, simultaneous=True), domain) != domain.zero
            for g in self.nonzero
        )


@dataclass(frozen=True)
class ConstructibleSet:
    """A finite union of locally closed subsets of one affine ambient space."""

    pieces: tuple[LocallyClosedSet, ...]

    def __post_init__(self):
        if any(not isinstance(piece, LocallyClosedSet) for piece in self.pieces):
            raise TypeError("all pieces must be locally closed sets")
        if self.pieces and any(
            piece.ideal.ring != self.pieces[0].ideal.ring for piece in self.pieces
        ):
            raise ValueError("constructible pieces must share an ambient polynomial ring")

    @property
    def ring(self):
        """The common polynomial ring, if this set has a nonempty representation."""
        return self.pieces[0].ideal.ring if self.pieces else None

    @property
    def closure_ideal(self):
        """Compute the Zariski closure by saturation and ideal intersection."""
        if not self.pieces:
            raise ValueError("an empty representation needs an ambient ring for its closure")
        result = self.ring.ideal(1)
        for piece in self.pieces:
            result = result.intersection(piece.closure_ideal)
        return result

    def contains(self, values, *, field=None):
        """Check membership in at least one of the locally closed pieces."""
        return any(piece.contains(values, field=field) for piece in self.pieces)

    def __or__(self, other):
        if not isinstance(other, ConstructibleSet) or (
            self.ring is not None and other.ring is not None and self.ring != other.ring
        ):
            return NotImplemented
        return ConstructibleSet(self.pieces + other.pieces)

    def __and__(self, other):
        if not isinstance(other, ConstructibleSet) or (
            self.ring is not None and other.ring is not None and self.ring != other.ring
        ):
            return NotImplemented
        return ConstructibleSet(
            tuple(
                LocallyClosedSet(a.ideal + b.ideal, a.nonzero + b.nonzero)
                for a in self.pieces
                for b in other.pieces
            )
        )

    def __sub__(self, other):
        """Compute set difference by distributing complements of locally closed pieces."""
        if not isinstance(other, ConstructibleSet) or (
            self.ring is not None and other.ring is not None and self.ring != other.ring
        ):
            return NotImplemented
        result = self
        for excluded in other.pieces:
            # Complement of V(I) intersect D(h) is union D(g_i) and V(h_j).
            complement = ConstructibleSet(
                tuple(
                    [
                        LocallyClosedSet(excluded.ideal.ring.ideal(), (g,))
                        for g in excluded.ideal.generators
                    ]
                    + [LocallyClosedSet(excluded.ideal.ring.ideal(h)) for h in excluded.nonzero]
                )
            )
            result = result & complement
        return result


@dataclass(frozen=True)
class ConstructibleImage:
    """An exact existential image of a polynomial or rational affine map.

    Membership is decided by solving the exact fiber equations over an
    algebraic closure. The representation is not yet a finite decomposition
    into locally closed subsets; ``closure_ideal`` is merely its closure.
    """

    mapping: object

    @property
    def closure_ideal(self):
        """Return the Zariski closure ideal, not a substitute for the exact image."""
        from .maps import RationalMap

        closure = self.mapping.image_closure
        ideal = closure if isinstance(self.mapping, RationalMap) else closure.ideal
        return ideal.radical()

    def as_constructible(self, *, max_minors=4096, max_exceptional_degree=4, max_branches=4096):
        """Return explicit exact locally closed image pieces over QQ.

        Uses structural linearization and rank conditions when possible,
        followed by constructive Chevalley elimination for nonlinear
        polynomial or rational maps. Branch limits raise explicitly.
        """
        from .chevalley import project_constructible
        from .linear_strata import linear_fiber_stratification
        from .maps import RationalMap
        from .schemes import AffineMorphism

        if isinstance(self.mapping, RationalMap):
            mapping = self.mapping
            if mapping.source.domain != sp.QQ:
                raise NotImplementedError("constructive projection currently requires QQ")
            ring = mapping.graph_ideal.ring
            graph = ring.ideal(
                *mapping.source.ideal.generators,
                *(
                    target * denominator.as_expr() - numerator.as_expr()
                    for target, (numerator, denominator) in zip(
                        mapping.targets, mapping._parts, strict=True
                    )
                ),
            )
            return project_constructible(
                ConstructibleSet((LocallyClosedSet(graph, mapping.denominator_factors),)),
                mapping.source.variables,
                max_branches=max_branches,
            )
        if not isinstance(self.mapping, AffineMorphism):
            raise TypeError("an affine polynomial morphism or rational map is required")
        mapping = self.mapping
        if mapping.source.domain != mapping.target.domain or mapping.source.domain != sp.QQ:
            raise NotImplementedError("image strata currently require rational coefficients")
        if set(mapping.source.variables) & set(mapping.target.variables):
            return self._chevalley_polynomial_image(max_branches=max_branches)
        if len(mapping.target.variables) == 1 and mapping.target.ideal.is_zero:
            return self._line_image_strata(max_exceptional_degree=max_exceptional_degree)
        variables = list(mapping.source.variables)
        target = mapping.target.variables
        expressions = [
            *mapping.source.ideal.generators,
            *(
                coordinate - expression
                for coordinate, expression in zip(target, mapping.coordinates, strict=True)
            ),
        ]
        # Unit-leading linear equations eliminate an entire source variable
        # without introducing exceptional parameter strata.
        while True:
            substitution = None
            for index, expression in enumerate(expressions):
                for variable in variables:
                    coefficient = sp.diff(expression, variable)
                    if coefficient in (sp.S.One, -sp.S.One):
                        substitution = (
                            index,
                            variable,
                            -expression.subs(variable, 0) / coefficient,
                        )
                        break
                if substitution is not None:
                    break
            if substitution is None:
                break
            index, variable, replacement = substitution
            expressions = [
                sp.expand(expr.subs(variable, replacement))
                for pos, expr in enumerate(expressions)
                if pos != index
            ]
            variables.remove(variable)
        if variables:
            try:
                family = linear_fiber_stratification(
                    expressions, tuple(variables), target, max_minors=max_minors
                )
            except NotImplementedError:
                return self._chevalley_polynomial_image(max_branches=max_branches)
        else:
            # After eliminating all source variables, the image is closed.
            return ConstructibleSet(
                (
                    LocallyClosedSet(
                        mapping.target.ambient.ideal(*expressions),
                    ),
                )
            )
        return ConstructibleSet(
            tuple(
                component
                for piece in family.pieces
                if piece.dimension >= 0
                for component in piece.locus.pieces
            )
        )

    def _chevalley_polynomial_image(self, *, max_branches=4096):
        """Project a polynomial graph with independent target-coordinate symbols."""
        from .chevalley import project_constructible

        mapping = self.mapping
        targets = mapping.target.variables
        shadows = tuple(sp.Dummy(f"image_{index}") for index in range(len(targets)))
        graph_ring = type(mapping.source.ambient)(
            (*mapping.source.variables, *shadows), sp.QQ, "lex"
        )
        graph = graph_ring.ideal(
            *mapping.source.ideal.generators,
            *(
                coordinate - expression
                for coordinate, expression in zip(shadows, mapping.coordinates, strict=True)
            ),
        )
        projected = project_constructible(
            ConstructibleSet((LocallyClosedSet(graph),)),
            mapping.source.variables,
            max_branches=max_branches,
        )
        rename = dict(zip(shadows, targets, strict=True))
        return ConstructibleSet(
            tuple(
                LocallyClosedSet(
                    mapping.target.ambient.ideal(
                        *(equation.xreplace(rename) for equation in piece.ideal.generators)
                    ),
                    tuple(g.xreplace(rename) for g in piece.nonzero),
                )
                for piece in projected.pieces
            )
        )

    def _line_image_strata(self, *, max_exceptional_degree=4):
        """Construct a finite/cofinite image of a polynomial map to A^1.

        Block-elimination Groebner leading coefficients identify a finite set
        outside which the specialized fiber basis remains nonunit. Exceptional
        irreducible factors are checked at one exact algebraic representative;
        Galois invariance makes the answer uniform over their conjugates.
        """
        mapping = self.mapping
        target = mapping.target.variables[0]
        ring = mapping.target.ambient
        closure = self.closure_ideal
        if not closure.is_zero:
            return ConstructibleSet((LocallyClosedSet(closure),))
        if not mapping.source.variables:
            raise ArithmeticError("dense line image cannot have zero source coordinates")
        equations = [
            *mapping.source.ideal.generators,
            *(
                coordinate - expression
                for coordinate, expression in zip(
                    mapping.target.variables, mapping.coordinates, strict=True
                )
            ),
        ]
        order = (*mapping.source.variables, target)
        basis = sp.groebner(equations, *order, order="lex", domain=sp.QQ)
        exceptional = sp.S.One
        for generator in basis.polys:
            leading = sp.Poly(generator.as_expr(), *mapping.source.variables).LC()
            exceptional *= leading
        polynomial = sp.Poly(exceptional, target, domain=sp.QQ).sqf_part()
        if polynomial.degree() <= 0:
            return ConstructibleSet((LocallyClosedSet(ring.ideal()),))
        factors = sp.factor_list(polynomial)[1]
        pieces = [LocallyClosedSet(ring.ideal(), (polynomial.as_expr(),))]
        for factor, _ in factors:
            if factor.degree() > max_exceptional_degree:
                raise NotImplementedError(
                    "exceptional factor degree exceeds max_exceptional_degree"
                )
            root = sp.rootof(factor.as_expr(), 0)
            field = sp.QQ.algebraic_field(root) if not root.is_Rational else sp.QQ
            if self.contains((root,), field=field):
                pieces.append(LocallyClosedSet(ring.ideal(factor.as_expr())))
        return ConstructibleSet(tuple(pieces))

    def strata(self, **options):
        """Return a finite locally closed image decomposition when computable."""
        return self.as_constructible(**options)

    def contains(self, values, *, field=None):
        """Decide geometric-fiber nonemptiness for an exact target point."""
        from .maps import RationalMap
        from .schemes import AffineMorphism

        mapping = self.mapping
        values = tuple(values)
        if isinstance(mapping, AffineMorphism):
            point = mapping.target.point(values, field=field)
            return not mapping.fiber(point).is_empty
        if not isinstance(mapping, RationalMap):
            raise TypeError("an affine polynomial morphism or rational map is required")
        if len(values) != len(mapping.targets):
            raise ValueError("one value is required per target coordinate")
        ring = mapping.source.ambient
        if field is not None:
            ring = type(ring)(ring.variables, field, ring.order)
        equations = list(mapping.source.ideal.generators)
        for value, (numerator, denominator) in zip(values, mapping._parts, strict=True):
            equations.append(numerator.as_expr() - value * denominator.as_expr())
        fiber_ideal = ring.ideal(*equations)
        if mapping.denominator_factors:
            fiber_ideal = fiber_ideal.saturate(sp.prod(mapping.denominator_factors))
        return not fiber_ideal.is_unit


__all__ = ["LocallyClosedSet", "ConstructibleSet", "ConstructibleImage"]
