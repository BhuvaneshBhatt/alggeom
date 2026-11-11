"""Explicit localization transitions between normalized projective-curve charts."""

from dataclasses import dataclass

import sympy as sp

from .schemes import AffineScheme, ProjectiveScheme


@dataclass(frozen=True)
class NormalizationPatch:
    """A normalized affine chart and its finite birational morphism."""

    chart: object
    mapping: object
    closure: object

    @property
    def normalized_scheme(self):
        """The normal affine source of the finite chart map."""
        return self.mapping.source


@dataclass(frozen=True)
class NormalizationTransition:
    """A regular transition on D(guard), with coordinates in its localization."""

    source_index: int
    target_index: int
    guard: sp.Expr
    coordinates: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class ProjectiveNormalizationAtlas:
    """Finite normalized charts, with explicit regular transition coordinates."""

    original: ProjectiveScheme
    patches: tuple[NormalizationPatch, ...]

    def _patch(self, index):
        return next((p for p in self.patches if p.chart.index == index), None)

    def overlap_guard(self, source_index, target_index):
        """Return the pullback of the target's distinguished homogeneous coordinate."""
        source = self._patch(source_index)
        if source is None or self._patch(target_index) is None:
            raise ValueError("the chart index is empty or not represented")
        coordinate = source.chart.embedding_coordinates[target_index]
        return source.mapping.pullback(coordinate).as_expr()

    def transition(self, source_index, target_index, *, max_localization_power=32):
        """Construct normalized chart transitions explicitly in A[guard^-1].

        The normalization's inverse rational functions express generators
        of the target patch's integral closure in its original curve field.
        Transporting these fractions through projective chart coordinates
        produces functions regular on the overlap. Ideal membership gives an
        explicit numerator divided by a power of the overlap guard; failure
        to establish regularity raises instead of returning a rational map.
        """
        source, target = self._patch(source_index), self._patch(target_index)
        if source is None or target is None:
            raise ValueError("the chart index is empty or not represented")
        guard = self.overlap_guard(source_index, target_index)
        target_original = dict(
            zip(
                target.chart.variables,
                source.chart.transition_coordinates(target.chart),
                strict=True,
            )
        )
        source_normal = dict(zip(source.chart.variables, source.mapping.coordinates, strict=True))
        ring = source.normalized_scheme.coordinate_ring
        localized = ring.localize(guard)
        images = []
        for coordinate in target.closure.inverse_functions:
            expression = coordinate.as_expr().subs(target_original, simultaneous=True)
            expression = sp.cancel(expression.subs(source_normal, simultaneous=True))
            numerator, denominator = sp.together(expression).as_numer_denom()
            # Find a polynomial representative after inverting the overlap
            # guard, including cases where the displayed denominator vanishes
            # on a removable-pole locus in the chart.
            ideal = ring.ideal + ring.ambient.ideal(denominator)
            power = sp.S.One
            for _ in range(max_localization_power + 1):
                certificate = ideal.membership(numerator * power)
                if certificate is not None:
                    representative = certificate[-1]
                    regular = sp.cancel(representative / power)
                    localized(regular)
                    images.append(regular)
                    break
                power *= guard
            else:
                raise NotImplementedError(
                    "could not express a transition coordinate on the overlap localization"
                )
        transition = NormalizationTransition(source_index, target_index, guard, tuple(images))
        field = ring.fraction_field()
        substitutions = dict(zip(target.normalized_scheme.variables, images, strict=True))
        for equation in target.normalized_scheme.ideal.generators:
            if field(equation.subs(substitutions, simultaneous=True)) != 0:
                raise ArithmeticError("normalization transition violates target chart equations")
        original_coordinates = source.chart.transition_coordinates(target.chart)
        original_substitutions = dict(
            zip(source.chart.variables, source.mapping.coordinates, strict=True)
        )
        for target_coordinate, original_coordinate in zip(
            target.mapping.coordinates, original_coordinates, strict=True
        ):
            computed = target_coordinate.subs(substitutions, simultaneous=True)
            expected = original_coordinate.subs(original_substitutions, simultaneous=True)
            if field(computed - expected) != 0:
                raise ArithmeticError(
                    "normalization transition disagrees with projective coordinates"
                )
        return transition

    def verify_cocycle(self, *, max_localization_power=32):
        """Verify identity and composition of transition maps in chart function fields.

        Algebraic equality on each nonempty triple overlap is checked in the
        normalized coordinate ring's exact fraction field. Regularity of the
        individual transitions was checked on their principal-open overlaps.
        """
        indices = tuple(p.chart.index for p in self.patches)
        transitions = {
            (i, j): self.transition(i, j, max_localization_power=max_localization_power)
            for i in indices
            for j in indices
        }
        for i in indices:
            source = self._patch(i)
            field = source.normalized_scheme.coordinate_ring.fraction_field()
            direct = transitions[(i, i)]
            if any(
                field(a - b) != 0
                for a, b in zip(direct.coordinates, source.normalized_scheme.variables, strict=True)
            ):
                return False
            for j in indices:
                middle = self._patch(j)
                mapping = dict(
                    zip(
                        middle.normalized_scheme.variables,
                        transitions[(i, j)].coordinates,
                        strict=True,
                    )
                )
                for k in indices:
                    composed = tuple(
                        c.subs(mapping, simultaneous=True) for c in transitions[(j, k)].coordinates
                    )
                    direct = transitions[(i, k)].coordinates
                    if any(field(a - b) != 0 for a, b in zip(composed, direct, strict=True)):
                        return False
        return True


def projective_normalization_atlas(curve, **limits):
    """Normalize standard projective curve charts and retain integral coordinates.

    Projective normalization is represented as an affine atlas with finite
    birational maps; each transition is computed and verified as a regular
    function on the localized overlap. No global projective embedding is
    constructed by this representation.
    """
    from .normalization import integral_closure

    if not isinstance(curve, ProjectiveScheme):
        raise TypeError("a projective scheme is required")
    if curve.dimension != 1 or not curve.is_integral:
        raise ValueError("projective normalization requires an integral projective curve")
    patches = []
    for index in range(len(curve.variables)):
        chart = curve.affine_chart(index)
        if not chart.scheme.is_empty:
            closure = integral_closure(chart.scheme, **limits)
            mapping = AffineScheme(closure.ring).morphism(
                chart.scheme, tuple(image.as_expr() for image in closure.mapping.images)
            )
            patches.append(NormalizationPatch(chart, mapping, closure))
    return ProjectiveNormalizationAtlas(curve, tuple(patches))


__all__ = [
    "NormalizationPatch",
    "NormalizationTransition",
    "ProjectiveNormalizationAtlas",
    "projective_normalization_atlas",
]
