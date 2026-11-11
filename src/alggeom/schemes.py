"""Affine and projective schemes, exact points, charts and polynomial morphisms."""

from dataclasses import dataclass
from functools import cached_property
from itertools import combinations

import sympy as sp
from algroots.domains import field_element
from algroots.rational_univariate.construction import rur_from_quotient

from ._algebra import coordinates, exact_field, integer
from .graded import degree, is_homogeneous, projective_closure
from .ideals import CoordinateRing, Ideal, PolynomialRing
from .rings import RingHomomorphism, _coefficient_embedding


def _ideal(value):
    if isinstance(value, CoordinateRing):
        return value.ideal
    if isinstance(value, PolynomialRing):
        return value.ideal()
    if not isinstance(value, Ideal):
        raise TypeError("an ideal, coordinate ring or polynomial ring is required")
    return value


def _same_ambient(left, right):
    if left.variables != right.variables or left.domain != right.domain:
        raise ValueError("schemes must have the same ambient coordinates and coefficient field")


def _point_coordinates(scheme, values, field):
    field = scheme.domain if field is None else exact_field(field)
    _coefficient_embedding(scheme.domain, field)
    values = tuple(field.to_sympy(field_element(v, field)) for v in values)
    if len(values) != len(scheme.variables):
        raise ValueError("one value is required for each ambient coordinate")
    substitution = dict(zip(scheme.variables, values, strict=True))
    if any(
        field_element(g.subs(substitution, simultaneous=True), field) != field.zero
        for g in scheme.ideal.generators
    ):
        raise ValueError("the coordinates do not satisfy the scheme equations")
    return values, field


def _value(expression, variables, point):
    substitution = dict(zip(variables, point.coordinates, strict=True))
    return point.field.to_sympy(
        field_element(sp.sympify(expression).subs(substitution, simultaneous=True), point.field)
    )


@dataclass(frozen=True)
class AffineScheme:
    """The affine scheme Spec(S/I), retaining nilpotents and the field of definition."""

    ideal: Ideal

    def __post_init__(self):
        object.__setattr__(self, "ideal", _ideal(self.ideal))

    @property
    def ambient(self):
        """The ambient polynomial ring."""
        return self.ideal.ring

    @property
    def variables(self):
        """The ordered ambient coordinate symbols."""
        return self.ambient.variables

    @property
    def domain(self):
        """The exact coefficient field."""
        return self.ambient.domain

    @cached_property
    def coordinate_ring(self):
        """The coordinate ring S/I, including its Krull and vector-space dimensions."""
        return CoordinateRing(self.ideal)

    @property
    def dimension(self):
        """The Krull dimension; the empty scheme has dimension -1."""
        return self.ideal.krull_dimension

    @property
    def is_empty(self):
        """Scheme emptiness, which is different from having no rational or real points."""
        return self.ideal.is_unit

    @property
    def is_reduced(self):
        """The field-relative is reduced property of this scheme."""
        from .properties import is_reduced

        return is_reduced(self)

    @property
    def is_irreducible(self):
        """The field-relative is irreducible property of this scheme."""
        from .properties import is_irreducible

        return is_irreducible(self)

    @property
    def is_integral(self):
        """The field-relative is integral property of this scheme."""
        from .properties import is_integral

        return is_integral(self)

    @property
    def is_equidimensional(self):
        """The field-relative is equidimensional property of this scheme."""
        from .properties import is_equidimensional

        return is_equidimensional(self)

    @property
    def is_unmixed(self):
        """The field-relative is unmixed property of this scheme."""
        from .properties import is_unmixed

        return is_unmixed(self)

    @property
    def codimension(self):
        """The field-relative codimension property of this scheme."""
        from .properties import codimension

        return codimension(self)

    def tangent_cone(self, point):
        """Return tangent cone at an exact point of this scheme."""
        from .local import tangent_cone

        return tangent_cone(self, point)

    def local_dimension(self, point):
        """Return local dimension at an exact point of this scheme."""
        from .local import local_dimension

        return local_dimension(self, point)

    def embedding_dimension(self, point):
        """Return embedding dimension at an exact point of this scheme."""
        from .local import embedding_dimension

        return embedding_dimension(self, point)

    def multiplicity(self, point):
        """Return multiplicity at an exact point of this scheme."""
        from .local import multiplicity

        return multiplicity(self, point)

    def tangent_space(self, point):
        """Return a basis of tangent directions in an affine chart."""
        from .maps import tangent_space

        return tangent_space(self, point)

    @property
    def is_normal(self):
        """Decide normality; projective schemes are checked on their affine charts."""
        from .normalization import is_normal

        return is_normal(self)

    def is_smooth(self, point=None):
        """Decide smoothness globally or at a specified point."""
        from .properties import is_smooth

        return is_smooth(self, point)

    def singular_locus(self, *, reduced=False):
        """Return the singular subscheme, with automatic dimension checks."""
        from .maps import singular_locus

        return singular_locus(self, reduced=reduced)

    def point_search(self, height=None, *, max_candidates=1000000):
        """Enumerate finite-field points or QQ points within an explicit height bound."""
        from .points import point_search

        return point_search(self, height=height, max_candidates=max_candidates)

    def point(self, values, *, field=None):
        """Validate a point over the base field or a specified exact extension."""
        return AffinePoint(self, values, field)

    def subscheme(self, *equations):
        """Impose additional equations with their scheme multiplicities."""
        return AffineScheme(self.coordinate_ring.quotient(*equations))

    def intersection(self, other):
        """Return the scheme-theoretic intersection in the same affine space."""
        _same_ambient(self, other)
        return AffineScheme(self.ideal + other.ideal)

    def union(self, other):
        """Return the scheme-theoretic union in the same affine space."""
        _same_ambient(self, other)
        return AffineScheme(self.ideal.intersection(other.ideal))

    def base_change(self, field):
        """Extend the exact coefficient field through its declared embedding."""
        return AffineScheme(self.coordinate_ring.base_change(field))

    def morphism(self, target, images):
        """Construct a polynomial morphism to an affine or projective target."""
        cls = ProjectiveMorphism if isinstance(target, ProjectiveScheme) else AffineMorphism
        return cls(self, target, images)

    def normalize(self, **limits):
        """Return the normalized affine scheme with its reduced integral-closure ring."""
        from .normalization import normalize

        return AffineScheme(normalize(self, **limits))

    def normalization_map(self, **limits):
        """Return the finite birational morphism from the normalization to this scheme."""
        from .normalization import normalization_map

        return normalization_map(self, **limits)

    def local_ring(self, point):
        """Return the local ring after extending coefficients to an exact point."""
        from .rational_functions import local_ring

        return local_ring(self, point)

    def noether_normalization(self, **options):
        """Construct a finite polynomial subring and the corresponding affine projection."""
        from .noether import noether_normalization

        return noether_normalization(self, **options)

    def projective_closure(self, *, homogenizer=None):
        """Return the saturated projective closure with one appended coordinate."""
        return ProjectiveScheme(projective_closure(self.ideal, homogenizer=homogenizer))

    def product(self, other, *, variables=None):
        """Form the affine product, with disjoint coordinates in factor order."""
        if not isinstance(other, AffineScheme) or self.domain != other.domain:
            raise ValueError("affine products require the same coefficient field")
        if variables is None:
            variables = tuple(sp.Dummy("left_" + v.name) for v in self.variables) + tuple(
                sp.Dummy("right_" + v.name) for v in other.variables
            )
        variables = coordinates(variables)
        if len(variables) != len(self.variables) + len(other.variables):
            raise ValueError("one product coordinate is required per factor coordinate")
        left = dict(zip(self.variables, variables[: len(self.variables)], strict=True))
        right = dict(zip(other.variables, variables[len(self.variables) :], strict=True))
        ring = PolynomialRing(variables, self.domain)
        return AffineScheme(
            ring.ideal(
                *(g.xreplace(left) for g in self.ideal.generators),
                *(g.xreplace(right) for g in other.ideal.generators),
            )
        )

    def principal_open(self, polynomial, *, inverse=None):
        """Realize D(f) as Spec(S[u]/(I,u*f-1)), retaining an explicit inverse coordinate."""
        f = self.coordinate_ring.element(polynomial).as_expr()
        inverse = sp.Dummy("inverse") if inverse is None else inverse
        ring = PolynomialRing((*self.variables, inverse), self.domain, self.ambient.order)
        return AffineScheme(ring.ideal(*self.ideal.generators, inverse * f - 1))

    def rational_points(self, *, field=None, max_dimension=256, max_candidates=1000000):
        """Enumerate all points in a specified field for a zero-dimensional scheme.

        Linear factors of an exact rational-univariate representation select
        precisely the points in that field. Multiplicity is retained in the
        scheme, not by repeating points. Positive-dimensional point searches
        and general Diophantine solving are outside this method.
        """
        scheme = self if field is None else self.base_change(field)
        from algroots.finite_fields import is_finite_field

        if is_finite_field(scheme.domain):
            from .points import point_search

            return tuple(
                self.point(p.coordinates, field=scheme.domain)
                for p in point_search(scheme, max_candidates=max_candidates)
            )
        if scheme.is_empty:
            return ()
        if scheme.dimension != 0:
            raise ValueError("complete point enumeration requires a zero-dimensional scheme")
        if not scheme.variables:
            return (self.point((), field=scheme.domain),)
        algebra = scheme.coordinate_ring.finite_algebra(max_dimension=max_dimension)
        representation = rur_from_quotient(algebra, sp.Dummy("point_parameter"))
        points = []
        for factor, _ in representation.defining_polynomial.factor_list()[1]:
            if factor.degree() != 1:
                continue
            root = -factor.nth(0) / factor.nth(1)
            denominator = representation.coordinate_denominator.eval(root)
            values = tuple(
                sp.cancel(numerator.eval(root) / denominator)
                for numerator in representation.coordinate_numerators
            )
            points.append(self.point(values, field=scheme.domain))
        return tuple(points)

    def reduced(self, **limits):
        """Return the reduced closed subscheme defined by the radical ideal."""
        from .components import radical

        return AffineScheme(radical(self.ideal, **limits))

    def irreducible_components(self, **limits):
        """Return reduced irreducible components over the declared coefficient field."""
        from .components import minimal_primes

        return tuple(AffineScheme(p) for p in minimal_primes(self.ideal, **limits))


@dataclass(frozen=True, eq=False)
class AffinePoint:
    """An exact field-valued point of an affine scheme, with a specified embedding."""

    scheme: AffineScheme
    coordinates: tuple[sp.Expr, ...]
    field: object = None

    __hash__ = None

    def __post_init__(self):
        if not isinstance(self.scheme, AffineScheme):
            raise TypeError("an affine point requires an AffineScheme")
        values, field = _point_coordinates(self.scheme, self.coordinates, self.field)
        object.__setattr__(self, "coordinates", values)
        object.__setattr__(self, "field", field)

    def __eq__(self, other):
        return (
            isinstance(other, AffinePoint)
            and self.scheme == other.scheme
            and self.field == other.field
            and all(
                field_element(a - b, self.field) == self.field.zero
                for a, b in zip(self.coordinates, other.coordinates, strict=True)
            )
        )

    def evaluate(self, expression):
        """Evaluate a regular coordinate-ring element at this point."""
        expression = self.scheme.coordinate_ring.element(expression).as_expr()
        return _value(expression, self.scheme.variables, self)

    def change_field(self, field):
        """Embed this point's values into a specified larger exact field."""
        _coefficient_embedding(self.field, exact_field(field))
        return AffinePoint(self.scheme, self.coordinates, field)


@dataclass(frozen=True)
class ProjectiveScheme:
    """A closed subscheme of ordinary projective space, with saturated homogeneous ideal."""

    ideal: Ideal

    def __post_init__(self):
        ideal = _ideal(self.ideal)
        if not ideal.ring.variables or not is_homogeneous(ideal):
            raise ValueError(
                "projective schemes require nonempty coordinates and a homogeneous ideal"
            )
        object.__setattr__(self, "ideal", ideal.saturate(ideal.ring.ideal(*ideal.ring.variables)))

    @property
    def ambient(self):
        """The homogeneous ambient polynomial ring."""
        return self.ideal.ring

    @property
    def variables(self):
        """The ordered homogeneous coordinates."""
        return self.ambient.variables

    @property
    def domain(self):
        """The coefficient field."""
        return self.ambient.domain

    @property
    def homogeneous_coordinate_ring(self):
        """The graded ring S/I, whose dimension differs from dimension of Proj."""
        return CoordinateRing(self.ideal)

    @property
    def dimension(self):
        """The geometric dimension of Proj; -1 denotes the empty scheme."""
        return max(-1, self.ideal.krull_dimension - 1)

    @property
    def degree(self):
        """The projective degree from the standard-graded Hilbert series."""
        return degree(self.ideal)

    @property
    def is_empty(self):
        """Whether the projective scheme is empty after irrelevant saturation."""
        return self.ideal.is_unit

    @property
    def is_reduced(self):
        """The field-relative is reduced property of this scheme."""
        from .properties import is_reduced

        return is_reduced(self)

    @property
    def is_irreducible(self):
        """The field-relative is irreducible property of this scheme."""
        from .properties import is_irreducible

        return is_irreducible(self)

    @property
    def is_integral(self):
        """The field-relative is integral property of this scheme."""
        from .properties import is_integral

        return is_integral(self)

    @property
    def is_equidimensional(self):
        """The field-relative is equidimensional property of this scheme."""
        from .properties import is_equidimensional

        return is_equidimensional(self)

    @property
    def is_unmixed(self):
        """The field-relative is unmixed property of this scheme."""
        from .properties import is_unmixed

        return is_unmixed(self)

    @property
    def codimension(self):
        """The field-relative codimension property of this scheme."""
        from .properties import codimension

        return codimension(self)

    def tangent_cone(self, point):
        """Return tangent cone at an exact point of this scheme."""
        from .local import tangent_cone

        return tangent_cone(self, point)

    def local_dimension(self, point):
        """Return local dimension at an exact point of this scheme."""
        from .local import local_dimension

        return local_dimension(self, point)

    def embedding_dimension(self, point):
        """Return embedding dimension at an exact point of this scheme."""
        from .local import embedding_dimension

        return embedding_dimension(self, point)

    def multiplicity(self, point):
        """Return multiplicity at an exact point of this scheme."""
        from .local import multiplicity

        return multiplicity(self, point)

    def tangent_space(self, point):
        """Return a basis of tangent directions in an affine chart."""
        from .maps import tangent_space

        return tangent_space(self, point)

    @property
    def is_normal(self):
        """Decide normality; projective schemes are checked on their affine charts."""
        from .normalization import is_normal

        return is_normal(self)

    def is_smooth(self, point=None):
        """Decide smoothness globally or at a specified point."""
        from .properties import is_smooth

        return is_smooth(self, point)

    def singular_locus(self, *, reduced=False):
        """Return the singular subscheme, with automatic dimension checks."""
        from .maps import singular_locus

        return singular_locus(self, reduced=reduced)

    def point_search(self, height=None, *, max_candidates=1000000):
        """Enumerate finite-field points or QQ points within an explicit height bound."""
        from .points import point_search

        return point_search(self, height=height, max_candidates=max_candidates)

    def point(self, values, *, field=None):
        """Validate homogeneous coordinates and normalize their first nonzero value."""
        return ProjectivePoint(self, values, field)

    def affine_chart(self, index, *, variables=None):
        """Dehomogenize on the chart where the indexed coordinate is nonzero."""
        return AffineChart(self, index, variables)

    def base_change(self, field):
        """Extend the coefficient field and saturate the resulting homogeneous ideal."""
        ring = self.homogeneous_coordinate_ring.base_change(field)
        return ProjectiveScheme(ring.ideal)

    def subscheme(self, *equations):
        """Impose additional homogeneous equations and saturate at the irrelevant ideal."""
        return ProjectiveScheme(self.ideal + self.ambient.ideal(*equations))

    def intersection(self, other):
        """Return the scheme-theoretic projective intersection."""
        _same_ambient(self, other)
        return ProjectiveScheme(self.ideal + other.ideal)

    def union(self, other):
        """Return the scheme-theoretic projective union."""
        _same_ambient(self, other)
        return ProjectiveScheme(self.ideal.intersection(other.ideal))

    def morphism(self, target, forms):
        """Construct a morphism given by homogeneous forms of equal degree."""
        return ProjectiveMorphism(self, target, forms)

    def rational_points(self, *, field=None, max_dimension=256, max_candidates=1000000):
        """Enumerate field-valued points of a zero-dimensional projective scheme by charts."""
        from algroots.finite_fields import is_finite_field

        scheme = self if field is None else self.base_change(field)
        if is_finite_field(scheme.domain):
            from .points import point_search

            return tuple(
                self.point(p.coordinates, field=scheme.domain)
                for p in point_search(scheme, max_candidates=max_candidates)
            )
        if self.is_empty:
            return ()
        if self.dimension != 0:
            raise ValueError("complete point enumeration requires a zero-dimensional scheme")
        result = []
        for index in range(len(self.variables)):
            chart = self.affine_chart(index)
            for point in chart.scheme.rational_points(field=field, max_dimension=max_dimension):
                image = chart.to_projective(point)
                if image not in result:
                    result.append(image)
        return tuple(result)

    def reduced(self, **limits):
        """Return the reduced projective closed subscheme."""
        from .components import radical

        return ProjectiveScheme(radical(self.ideal, **limits))

    def irreducible_components(self, **limits):
        """Return the reduced projective components over the declared coefficient field."""
        from .components import minimal_primes

        return tuple(ProjectiveScheme(p) for p in minimal_primes(self.ideal, **limits))


@dataclass(frozen=True, eq=False)
class ProjectivePoint:
    """A projective point over an exact field, with coordinates canonical up to scaling."""

    scheme: ProjectiveScheme
    coordinates: tuple[sp.Expr, ...]
    field: object = None

    __hash__ = None

    def __post_init__(self):
        if not isinstance(self.scheme, ProjectiveScheme):
            raise TypeError("a projective point requires a ProjectiveScheme")
        values, field = _point_coordinates(self.scheme, self.coordinates, self.field)
        native = tuple(field_element(v, field) for v in values)
        nonzero = next((v for v in native if v != field.zero), None)
        if nonzero is None:
            raise ValueError("projective coordinates cannot all be zero")
        object.__setattr__(self, "coordinates", tuple(field.to_sympy(v / nonzero) for v in native))
        object.__setattr__(self, "field", field)

    def __eq__(self, other):
        return (
            isinstance(other, ProjectivePoint)
            and self.scheme == other.scheme
            and self.field == other.field
            and all(
                field_element(a - b, self.field) == self.field.zero
                for a, b in zip(self.coordinates, other.coordinates, strict=True)
            )
        )


@dataclass(frozen=True)
class AffineChart:
    """A projective affine patch with its embedding and transition coordinates."""

    projective_scheme: ProjectiveScheme
    index: int
    variables: tuple[sp.Symbol, ...] = None

    def __post_init__(self):
        if not isinstance(self.projective_scheme, ProjectiveScheme):
            raise TypeError("an affine chart requires a ProjectiveScheme")
        index = integer(self.index, "index")
        if index >= len(self.projective_scheme.variables):
            raise ValueError("the chart index is outside the homogeneous coordinates")
        remaining = tuple(v for i, v in enumerate(self.projective_scheme.variables) if i != index)
        variables = remaining if self.variables is None else coordinates(self.variables)
        if len(variables) != len(remaining):
            raise ValueError(
                "one chart coordinate is required per remaining homogeneous coordinate"
            )
        object.__setattr__(self, "variables", variables)

    @cached_property
    def embedding_coordinates(self):
        """Homogeneous coordinate expressions for the affine-chart embedding."""
        values = iter(self.variables)
        return tuple(
            sp.S.One if i == self.index else next(values)
            for i in range(len(self.projective_scheme.variables))
        )

    @cached_property
    def scheme(self):
        """The affine scheme obtained by setting the selected homogeneous coordinate to one."""
        substitution = dict(
            zip(self.projective_scheme.variables, self.embedding_coordinates, strict=True)
        )
        ring = PolynomialRing(self.variables, self.projective_scheme.domain)
        return AffineScheme(
            ring.ideal(
                *(
                    g.subs(substitution, simultaneous=True)
                    for g in self.projective_scheme.ideal.generators
                )
            )
        )

    @cached_property
    def embedding(self):
        """The polynomial morphism from this chart to its projective scheme."""
        return ProjectiveMorphism(self.scheme, self.projective_scheme, self.embedding_coordinates)

    def to_projective(self, point):
        """Embed a point of this affine chart into projective space."""
        return self.embedding(point)

    def from_projective(self, point):
        """Convert a point in this open chart to its affine coordinates."""
        if not isinstance(point, ProjectivePoint) or point.scheme != self.projective_scheme:
            raise ValueError("the point must belong to the chart's projective scheme")
        denominator = field_element(point.coordinates[self.index], point.field)
        if denominator == point.field.zero:
            raise ValueError("the point does not belong to this affine chart")
        values = tuple(
            point.field.to_sympy(field_element(v, point.field) / denominator)
            for i, v in enumerate(point.coordinates)
            if i != self.index
        )
        return self.scheme.point(values, field=point.field)

    def transition_coordinates(self, other):
        """Return rational coordinates into another chart on their overlap."""
        if not isinstance(other, AffineChart) or other.projective_scheme != self.projective_scheme:
            raise ValueError("chart transitions require the same projective scheme")
        denominator = self.embedding_coordinates[other.index]
        return tuple(
            sp.cancel(v / denominator)
            for i, v in enumerate(self.embedding_coordinates)
            if i != other.index
        )

    def overlap_guard(self, other):
        """The nonvanishing condition defining overlap with another chart."""
        self.transition_coordinates(other)
        return sp.Ne(self.embedding_coordinates[other.index], 0)


@dataclass(frozen=True)
class AffineMorphism:
    """A polynomial morphism of affine schemes over one coefficient field."""

    source: AffineScheme
    target: AffineScheme
    coordinates: tuple[sp.Expr, ...]

    def __post_init__(self):
        if not isinstance(self.source, AffineScheme) or not isinstance(self.target, AffineScheme):
            raise TypeError("an affine morphism requires affine source and target schemes")
        images = tuple(self.source.coordinate_ring.element(v).as_expr() for v in self.coordinates)
        object.__setattr__(self, "coordinates", images)
        object.__setattr__(
            self,
            "_pullback",
            RingHomomorphism(self.target.coordinate_ring, self.source.coordinate_ring, images),
        )

    @property
    def pullback(self):
        """The contravariant map from target to source coordinate rings."""
        return self._pullback

    def __call__(self, point):
        if not isinstance(point, AffinePoint) or point.scheme != self.source:
            raise ValueError("the point must belong to the source scheme")
        return self.target.point(
            tuple(_value(f, self.source.variables, point) for f in self.coordinates),
            field=point.field,
        )

    def compose(self, inner):
        """Return self after inner as a polynomial morphism."""
        if not isinstance(inner, AffineMorphism) or inner.target != self.source:
            raise ValueError("the intermediate schemes must agree")
        return AffineMorphism(
            inner.source, self.target, tuple(inner.pullback(f).as_expr() for f in self.coordinates)
        )

    @property
    def graph_ideal(self):
        """The graph ideal, with fresh target symbols following the source coordinates."""
        return self.pullback.graph_ideal

    @property
    def image_closure(self):
        """The scheme-theoretic image closure, computed from the pullback kernel."""
        return AffineScheme(self.pullback.kernel)

    @property
    def is_dominant(self):
        """Decide whether the image is topologically dense in the target."""
        from .morphism_invariants import is_dominant

        return is_dominant(self)

    @property
    def is_finite(self):
        """Decide module finiteness over the affine target coordinate ring."""
        from .morphism_invariants import is_finite

        return is_finite(self)

    @property
    def generic_degree(self):
        """Return generic fiber length for a dominant generically finite affine morphism."""
        from .morphism_invariants import morphism_degree

        return morphism_degree(self)

    def ramification_locus(self):
        """Return the non-unramified subscheme of a finite affine morphism."""
        from .morphism_invariants import ramification_locus

        return ramification_locus(self)

    def branch_locus(self):
        """Return the image closure of the relative non-unramified locus."""
        from .morphism_invariants import branch_locus

        return branch_locus(self)

    @property
    def constructible_image(self):
        """Return the exact geometric image with an existential fiber predicate."""
        from .constructible import ConstructibleImage

        return ConstructibleImage(self)

    def preimage(self, subscheme):
        """Return the scheme-theoretic inverse image of a target closed subscheme."""
        if not isinstance(subscheme, AffineScheme):
            raise TypeError("an affine closed subscheme is required")
        _same_ambient(self.target, subscheme)
        if not self.target.ideal.is_subset(subscheme.ideal):
            raise ValueError("the subscheme is not contained in the target")
        return self.source.subscheme(*(self.pullback(g) for g in subscheme.ideal.generators))

    def fiber(self, point):
        """Return the exact fiber over a field-valued target point, including nilpotents."""
        if not isinstance(point, AffinePoint) or point.scheme != self.target:
            raise ValueError("the point must belong to the target scheme")
        source = self.source.base_change(point.field)
        return source.subscheme(
            *(f - v for f, v in zip(self.coordinates, point.coordinates, strict=True))
        )

    def fiber_product(self, other, *, variables=None):
        """Form the affine fiber product of two morphisms to the same target."""
        if not isinstance(other, AffineMorphism) or other.target != self.target:
            raise ValueError("fiber products require the same target scheme")
        product = self.source.product(other.source, variables=variables)
        left = dict(
            zip(self.source.variables, product.variables[: len(self.source.variables)], strict=True)
        )
        right = dict(
            zip(
                other.source.variables, product.variables[len(self.source.variables) :], strict=True
            )
        )
        return product.subscheme(
            *(
                f.xreplace(left) - g.xreplace(right)
                for f, g in zip(self.coordinates, other.coordinates, strict=True)
            )
        )

    def base_change(self, field):
        """Extend source, target and coordinate expressions to a larger exact field."""
        return AffineMorphism(
            self.source.base_change(field), self.target.base_change(field), self.coordinates
        )

    def inverse(self):
        """Return the inverse when the coordinate-ring pullback is an isomorphism."""
        inverse = self.pullback.inverse()
        return AffineMorphism(self.target, self.source, tuple(v.as_expr() for v in inverse.images))


@dataclass(frozen=True)
class ProjectiveMorphism:
    """A morphism to ordinary projective space represented by global polynomial coordinates.

    A projective source requires homogeneous forms of equal degree. Coordinates
    must have no common zero on the source. Affine sources allow arbitrary
    polynomial coordinates with no common affine zero. This presentation does
    not represent every line bundle or every abstract scheme morphism.
    """

    source: AffineScheme | ProjectiveScheme
    target: ProjectiveScheme
    coordinates: tuple[sp.Expr, ...]

    def __post_init__(self):
        if not isinstance(self.source, (AffineScheme, ProjectiveScheme)) or not isinstance(
            self.target, ProjectiveScheme
        ):
            raise TypeError("a projective target and an affine or projective source are required")
        if self.source.domain != self.target.domain:
            raise ValueError("source and target must have the same coefficient field")
        forms = tuple(self.source.ambient.convert(f) for f in self.coordinates)
        if len(forms) != len(self.target.variables):
            raise ValueError("one coordinate form is required per target coordinate")
        if isinstance(self.source, ProjectiveScheme):
            degrees = [
                {sum(m) for m in self.source.ambient._poly(f).monoms()} for f in forms if f != 0
            ]
            if any(len(d) != 1 for d in degrees) or len(set().union(*degrees)) > 1:
                raise ValueError(
                    "projective coordinate forms must be homogeneous of the same degree"
                )
        object.__setattr__(self, "coordinates", forms)
        substitution = dict(zip(self.target.variables, forms, strict=True))
        if any(
            not self.source.ideal.contains(g.subs(substitution, simultaneous=True))
            for g in self.target.ideal.generators
        ):
            raise ValueError("the coordinate forms do not satisfy the target equations")
        base = self.source.ideal + self.source.ambient.ideal(*forms)
        if isinstance(self.source, ProjectiveScheme):
            base = base.saturate(self.source.ambient.ideal(*self.source.variables))
        if not base.is_unit:
            raise ValueError("the coordinate forms have a base locus on the source")

    def __call__(self, point):
        if not isinstance(point, (AffinePoint, ProjectivePoint)) or point.scheme != self.source:
            raise ValueError("the point must belong to the source scheme")
        return self.target.point(
            tuple(_value(f, self.source.variables, point) for f in self.coordinates),
            field=point.field,
        )

    def compose(self, inner):
        """Return self after another morphism into its projective source."""
        if not isinstance(inner, ProjectiveMorphism) or inner.target != self.source:
            raise ValueError("the intermediate projective schemes must agree")
        substitution = dict(zip(self.source.variables, inner.coordinates, strict=True))
        return ProjectiveMorphism(
            inner.source,
            self.target,
            tuple(f.subs(substitution, simultaneous=True) for f in self.coordinates),
        )

    def fiber(self, point):
        """Return the fiber over a projective point, as an affine or projective scheme."""
        if not isinstance(point, ProjectivePoint) or point.scheme != self.target:
            raise ValueError("the point must belong to the target scheme")
        source = self.source.base_change(point.field)
        pivot = next(i for i, v in enumerate(point.coordinates) if v != 0)
        return source.subscheme(
            *(
                point.coordinates[pivot] * f - v * self.coordinates[pivot]
                for f, v in zip(self.coordinates, point.coordinates, strict=True)
            )
        )

    def preimage(self, subscheme):
        """Pull back the homogeneous equations of a projective closed subscheme."""
        if not isinstance(subscheme, ProjectiveScheme):
            raise TypeError("a projective closed subscheme is required")
        _same_ambient(self.target, subscheme)
        if not self.target.ideal.is_subset(subscheme.ideal):
            raise ValueError("the subscheme is not contained in the target")
        substitution = dict(zip(self.target.variables, self.coordinates, strict=True))
        return self.source.subscheme(
            *(g.subs(substitution, simultaneous=True) for g in subscheme.ideal.generators)
        )

    @cached_property
    def graph_ideal(self):
        """The graph ideal in the product coordinates, saturated in both projective blocks."""
        target = tuple(sp.Dummy("target_" + v.name) for v in self.target.variables)
        ring = PolynomialRing((*self.source.variables, *target), self.source.domain, "lex")
        substitution = dict(zip(self.target.variables, target, strict=True))
        equations = [
            *self.source.ideal.generators,
            *(g.xreplace(substitution) for g in self.target.ideal.generators),
        ]
        equations.extend(
            target[i] * self.coordinates[j] - target[j] * self.coordinates[i]
            for i, j in combinations(range(len(target)), 2)
        )
        ideal = ring.ideal(*equations).saturate(ring.ideal(*target))
        if isinstance(self.source, ProjectiveScheme):
            ideal = ideal.saturate(ring.ideal(*self.source.variables))
        return ideal

    @property
    def image_closure(self):
        """The projective scheme-theoretic image closure by saturated graph elimination."""
        ideal = self.graph_ideal.eliminate(*self.source.variables)
        rename = dict(zip(ideal.ring.variables, self.target.variables, strict=True))
        return ProjectiveScheme(
            self.target.ambient.ideal(*(g.xreplace(rename) for g in ideal.generators))
        )

    def base_change(self, field):
        """Extend the coefficient field of this polynomial presentation."""
        return ProjectiveMorphism(
            self.source.base_change(field), self.target.base_change(field), self.coordinates
        )


__all__ = [
    "AffineScheme",
    "ProjectiveScheme",
    "AffinePoint",
    "ProjectivePoint",
    "AffineChart",
    "AffineMorphism",
    "ProjectiveMorphism",
]
