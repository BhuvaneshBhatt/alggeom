"""Graphs, image closures and affine geometric diagnostics for exact polynomial ideals."""

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations

import sympy as sp
from algroots.domains import domain_matrix, field_element

from ._algebra import coordinates, rational_parts, visible_denominators
from .ideals import CoordinateRing, Ideal, PolynomialRing


def _point(ideal, point):
    variables, domain = ideal.ring.variables, ideal.ring.domain
    if isinstance(point, Mapping):
        if set(point) != set(variables):
            raise ValueError("a point must assign every coordinate exactly once")
        values = tuple(point[v] for v in variables)
    else:
        values = tuple(point)
        if len(values) != len(variables):
            raise ValueError("one point coordinate is required per ambient variable")
    values = tuple(domain.to_sympy(field_element(v, domain)) for v in values)
    substitution = dict(zip(variables, values, strict=True))
    if any(
        field_element(g.subs(substitution, simultaneous=True), domain) != domain.zero
        for g in ideal.generators
    ):
        raise ValueError("the supplied point does not lie on the source scheme")
    return substitution


@dataclass(frozen=True)
class RationalMap:
    """Rational expressions on an explicit principal open subset of an affine scheme.

    Original denominator exclusions are retained even after cancellation.
    Graph and image ideals describe Zariski closures. They do not encode an
    exact constructible image, a real locus, or selected real branches.
    """

    source: CoordinateRing
    coordinates: tuple[sp.Expr, ...]
    targets: tuple[sp.Symbol, ...]
    excluded: tuple[sp.Expr, ...] = ()

    def __post_init__(self):
        source = self.source
        if isinstance(source, PolynomialRing):
            source = CoordinateRing(source.ideal())
        elif isinstance(source, Ideal):
            source = CoordinateRing(source)
        if not isinstance(source, CoordinateRing):
            raise TypeError("source must be a polynomial ring, ideal, or coordinate ring")
        object.__setattr__(self, "source", source)
        object.__setattr__(
            self,
            "coordinates",
            tuple(
                sp.sympify(f.as_expr() if isinstance(f, sp.Poly) else f) for f in self.coordinates
            ),
        )
        object.__setattr__(self, "targets", coordinates(self.targets))
        object.__setattr__(
            self, "excluded", tuple(source.ambient.convert(f) for f in self.excluded)
        )
        parameters = source.domain.symbols if source.domain.is_FractionField else ()
        if len(self.coordinates) != len(self.targets) or set(self.targets) & {
            *source.variables,
            *parameters,
        }:
            raise ValueError("one fresh target Symbol is required per map coordinate")
        for f in self.coordinates:
            rational_parts(f, source.ambient)

    @cached_property
    def _parts(self):
        return tuple(rational_parts(f, self.source.ambient) for f in self.coordinates)

    @property
    def functions(self):
        """Reduced coordinate expressions over the ambient polynomial fraction field."""
        return tuple(p.as_expr() / q.as_expr() for p, q in self._parts)

    @cached_property
    def denominator_factors(self):
        """Polynomial nonvanishing factors defining the declared affine input domain."""
        factors = (*self.excluded, *(d for f in self.coordinates for d in visible_denominators(f)))
        return tuple(dict.fromkeys(self.source.ambient.convert(f) for f in factors))

    @property
    def domain_guard(self):
        """Nonvanishing conditions, independent of the source's defining equations."""
        return sp.And(*(sp.Ne(f, 0) for f in self.denominator_factors))

    @cached_property
    def graph_ideal(self):
        """The closure of the graph over its declared denominator-free input domain."""
        ring = PolynomialRing((*self.source.variables, *self.targets), self.source.domain, "lex")
        equations = (
            *self.source.ideal.generators,
            *(
                v * q.as_expr() - p.as_expr()
                for v, (p, q) in zip(self.targets, self._parts, strict=True)
            ),
        )
        ideal = ring.ideal(*equations)
        denominator = sp.prod(self.denominator_factors)
        return ideal.saturate(denominator)

    @cached_property
    def image_closure(self):
        """The defining ideal of the scheme-theoretic image closure in target space."""
        return self.graph_ideal.eliminate(*self.source.variables)

    @property
    def constructible_image(self):
        """Return the exact geometric image, retaining excluded denominator loci."""
        from .constructible import ConstructibleImage

        return ConstructibleImage(self)

    def evaluate(self, point):
        """Evaluate an exact coefficient-field point, checking equations and input poles."""
        substitution = _point(self.source.ideal, point)
        domain = self.source.domain
        if any(
            field_element(f.subs(substitution, simultaneous=True), domain) == domain.zero
            for f in self.denominator_factors
        ):
            raise ValueError("the point is excluded by an input denominator")
        return tuple(
            domain.to_sympy(field_element(f.subs(substitution, simultaneous=True), domain))
            for f in self.functions
        )

    def compose(self, inner, *, targets=None):
        """Compose maps while retaining poles from both stages.

        The inner image must lie in the outer source. New output symbols are
        chosen when the outer targets would capture an inner source coordinate.
        """
        if (
            not isinstance(inner, RationalMap)
            or len(self.source.variables) != len(inner.targets)
            or self.source.domain != inner.source.domain
        ):
            raise ValueError(
                "composition requires matching coefficient fields and source/target arities"
            )
        correspondence = dict(zip(self.source.variables, inner.targets, strict=True))
        if any(
            not inner.graph_ideal.contains(g.xreplace(correspondence))
            for g in self.source.ideal.generators
        ):
            raise ValueError("the inner image does not lie in the outer source")
        substitution = dict(zip(self.source.variables, inner.functions, strict=True))
        exclusions = list(inner.denominator_factors)
        # Check pulled-back poles before evaluating the outer functions. A constant
        # map into a pole has an empty composite domain, not a value of zoo.
        for factor in self.denominator_factors:
            pulled_back = factor.subs(substitution, simultaneous=True)
            numerator, denominator = rational_parts(pulled_back, inner.source.ambient)
            if inner.source.ideal.contains(numerator.as_expr()):
                raise ValueError("composite has an empty domain: outer pole on inner image")
            exclusions.extend((numerator.as_expr(), denominator.as_expr()))
        expressions = tuple(f.subs(substitution, simultaneous=True) for f in self.functions)
        if targets is None:
            targets = tuple(
                sp.Dummy(v.name) if v in inner.source.variables else v for v in self.targets
            )
        return RationalMap(inner.source, expressions, targets, tuple(exclusions))


def _map(expressions, parameters, targets, source, domain):
    parameters = coordinates(parameters)
    if source is None:
        source = PolynomialRing(parameters, sp.QQ if domain is None else domain)
    elif isinstance(source, Ideal):
        source = CoordinateRing(source)
    if not isinstance(source, (PolynomialRing, CoordinateRing)):
        raise TypeError("source must be a polynomial ring, ideal, or coordinate ring")
    if source.variables != parameters or (
        domain is not None and source.domain != PolynomialRing(parameters, domain).domain
    ):
        raise ValueError(
            "source coordinates and coefficient field must match the declared parameters"
        )
    expressions = tuple(expressions)
    targets = (
        tuple(sp.Dummy(f"image_{i}") for i in range(len(expressions)))
        if targets is None
        else targets
    )
    return RationalMap(source, expressions, targets)


def implicitize(expressions, parameters, *, targets=None, source=None, domain=None):
    """Return the image-closure ideal, inferring the field from source or defaulting to QQ."""
    return _map(expressions, parameters, targets, source, domain).image_closure


def graph_ideal(expressions, parameters, *, targets=None, source=None, domain=None):
    """Return the saturated graph closure, with source coordinates preceding targets.

    The field is inferred from source, or defaults to QQ when no source is given.
    """
    return _map(expressions, parameters, targets, source, domain).graph_ideal


def jacobian(value):
    """Return the defining-equation Jacobian of an ideal, ring or embedded scheme."""
    from .properties import defining_ideal

    ideal = defining_ideal(value)
    return sp.ImmutableMatrix(
        len(ideal.generators),
        len(ideal.ring.variables),
        lambda i, j: sp.diff(ideal.generators[i], ideal.ring.variables[j]),
    )


def tangent_space(value, point):
    """Return a column basis of tangent directions after extending the point field.

    Projective inputs use the affine chart of the first nonzero point coordinate.
    Nonreduced and mixed-dimensional schemes retain their defining equations.
    """
    from .local import _affine_point

    ideal, values = _affine_point(value, point)
    matrix = jacobian(ideal).subs(
        dict(zip(ideal.ring.variables, values, strict=True)), simultaneous=True
    )
    return sp.ImmutableMatrix(
        domain_matrix(matrix, ideal.ring.domain).nullspace(divide_last=True).to_Matrix().T
    )


def tangent_ideal(value, point):
    """Return the embedded affine tangent space, translated through the chosen point."""
    from .local import _affine_point

    ideal, values = _affine_point(value, point)
    matrix = jacobian(ideal).subs(
        dict(zip(ideal.ring.variables, values, strict=True)), simultaneous=True
    )
    displacement = sp.Matrix([v - p for v, p in zip(ideal.ring.variables, values, strict=True)])
    return ideal.ring.ideal(*list(matrix * displacement))


def _jacobian_ideal(ideal, codimension):
    matrix = jacobian(ideal)
    return ideal + ideal.ring.ideal(
        *(
            matrix.extract(rows, columns).det()
            for rows in combinations(range(matrix.rows), codimension)
            for columns in combinations(range(matrix.cols), codimension)
        )
    )


def singular_locus(value, *, reduced=False, assume_equidimensional=False):
    """Compute a singular subscheme, inferring equidimensionality when necessary.

    Equidimensional inputs retain the Jacobian scheme structure unless reduced
    is requested. Mixed-dimensional inputs require reduced=True: component
    dimension strata are saturated away from higher-dimensional components,
    then their singular sets are joined. Scheme inputs return schemes;
    projective results are saturated at the irrelevant ideal.
    """
    from .properties import defining_ideal
    from .schemes import AffineScheme, ProjectiveScheme

    ideal = defining_ideal(value)
    if not isinstance(reduced, bool) or not isinstance(assume_equidimensional, bool):
        raise ValueError("reduced and assume_equidimensional must be Boolean")
    codimension = len(ideal.ring.variables) - ideal.krull_dimension
    if ideal.is_unit:
        result = ideal
    elif (
        ideal.krull_dimension == 0 or len(ideal.generators) == codimension or assume_equidimensional
    ):
        result = _jacobian_ideal(ideal, codimension)
    else:
        primes = ideal.minimal_primes()
        dimensions = {p.krull_dimension for p in primes}
        if len(dimensions) <= 1:
            result = _jacobian_ideal(ideal, codimension)
        else:
            if not reduced:
                raise ValueError("mixed-dimensional singular loci require reduced=True")
            result = ideal.ring.ideal(1)
            for dimension in sorted(dimensions):
                current = ideal.ring.ideal(1)
                higher = ideal.ring.ideal(1)
                for prime in primes:
                    if prime.krull_dimension == dimension:
                        current = current.intersection(prime)
                    elif prime.krull_dimension > dimension:
                        higher = higher.intersection(prime)
                locus = (
                    current + _jacobian_ideal(ideal, len(ideal.ring.variables) - dimension)
                ).saturate(higher)
                result = result.intersection(locus)
    if reduced:
        result = result.radical()
    if isinstance(value, ProjectiveScheme):
        return ProjectiveScheme(result)
    return AffineScheme(result) if isinstance(value, AffineScheme) else result


__all__ = [
    "RationalMap",
    "implicitize",
    "graph_ideal",
    "jacobian",
    "tangent_space",
    "tangent_ideal",
    "singular_locus",
]
