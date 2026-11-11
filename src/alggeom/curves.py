"""Plane-curve genus, normalization analysis and genus-zero parametrization.

Models are reduced irreducible plane equations over QQ or number fields.
Normalization analysis checks geometric integrality and resolves singularities
to compute the smooth projective genus and adjoints. It does not construct a
general affine integral closure. Rational curves additionally have explicit
birational parametrizations and projective normalization morphisms.
"""

from dataclasses import dataclass
from functools import cached_property

import sympy as sp
from algroots.curves import PlaneCurve as _PlaneModel

from .ideals import PolynomialRing
from .maps import RationalMap
from .schemes import AffineScheme, ProjectiveScheme


@dataclass(frozen=True)
class PlaneCurve:
    """A reduced irreducible affine plane model over QQ or a number field.

    The affine equation and its projective closure share a function field.
    Genus refers to the smooth projective normalization, so singular plane
    degree alone does not determine it. Geometrically reducible models are
    rejected by normalization analysis. Rational-function constant fields,
    arbitrary space curves and higher-dimensional normalization are outside
    this interface's scope.
    """

    equation: sp.Expr
    variables: tuple[sp.Symbol, sp.Symbol]
    domain: object = sp.QQ

    def __post_init__(self):
        model = _PlaneModel.from_polynomial(self.equation, self.variables, domain=self.domain)
        object.__setattr__(self, "equation", model.polynomial.as_expr())
        object.__setattr__(self, "variables", model.variables)
        object.__setattr__(self, "domain", model.domain)
        object.__setattr__(self, "_model", model)

    @property
    def degree(self):
        """The total degree of the plane model, not the normalized curve's genus."""
        return self._model.degree

    @cached_property
    def affine_scheme(self):
        """The affine closed scheme defined by this equation in A^2."""
        return AffineScheme(PolynomialRing(self.variables, self.domain).ideal(self.equation))

    @cached_property
    def projective_scheme(self):
        """The projective closure, with a fresh final homogenizing coordinate."""
        return self.affine_scheme.projective_closure()

    @property
    def genus(self):
        """The genus of the smooth projective normalization over the declared constants."""
        return self.normalization().genus

    def normalization(self, *, max_steps=256, max_dimension=256):
        """Return algroots' exact CurveNormalization analysis certificate.

        It records genus, arithmetic_genus, delta_sum, singularities including
        infinitely near points, and canonical_adjoints, with verify() replay.
        This is resolution and adjoint data, not an integral-closure coordinate
        ring. Computation limits raise without a partial genus claim.
        """
        return self._model.normalization(max_steps=max_steps, max_dimension=max_dimension)

    def rationality(
        self,
        parameter=None,
        *,
        point=None,
        extend_constants=False,
        max_steps=256,
        max_dimension=256,
    ):
        """Return an exact CurveRationalityResult with maps or a rationality obstruction.

        Status is rational, rational_over_extension, positive_genus, or
        no_ground_field_point. A genus-zero obstruction over QQ is decided by
        a conic model. Other number fields need a supplied point witness or
        extend_constants=True; failure to find a point is not called an
        obstruction. Any constants extension and the inverse birational map
        remain available in this specialist certificate.
        """
        return self._model.rational_parametrization(
            parameter,
            point=point,
            extend_constants=extend_constants,
            max_steps=max_steps,
            max_dimension=max_dimension,
        )

    def parametrize(self, parameter=None, **options):
        """Return a RationalMap from an open subset of A^1 onto this plane curve.

        Denominator exclusions describe this affine parameter chart. Use
        rationality() for inverse coordinates and exceptional points, or
        normalization_map() to include the projective parameter at infinity.
        A positive genus or missing ground-field point raises ValueError.
        """
        result = self.rationality(parameter, **options)
        if result.parametrization is None:
            raise ValueError(f"the curve has no requested parametrization: {result.status}")
        data = result.parametrization
        source = PolynomialRing((data.parameter,), data.domain)
        return RationalMap(source, data.affine_coordinates, self.variables)

    def affine_normalization(self, **limits):
        """Return the affine integral-closure scheme, including positive-genus plane models."""
        return self.affine_scheme.normalize(**limits)

    def projective_normalization_atlas(self, **limits):
        """Return exact affine normalization maps covering the projective model.

        The returned atlas records normalized affine patches and overlap opens.
        It does not provide explicit transition coordinates for newly adjoined
        integral generators or a single projective source embedding.
        """
        from .projective_normalization import projective_normalization_atlas

        return projective_normalization_atlas(self.projective_scheme, **limits)

    def normalization_map(self, *, parameters=None, **options):
        """Return the projective normalization morphism P^1 -> C for a rational curve.

        The target is base changed if extend_constants=True supplies a
        parametrization over a larger field. Positive-genus curves currently
        supply normalization analysis only, without a normalized scheme model.
        """
        result = self.rationality(**options)
        if result.parametrization is None:
            raise ValueError(f"the curve has no requested parametrization: {result.status}")
        data = result.parametrization
        parameters = (sp.Dummy("u"), sp.Dummy("v")) if parameters is None else tuple(parameters)
        if len(parameters) != 2:
            raise ValueError("two homogeneous parameter coordinates are required")
        source = ProjectiveScheme(PolynomialRing(parameters, data.domain).ideal())
        target = self.projective_scheme.base_change(data.domain)
        return source.morphism(target, data.homogeneous_coordinates(*parameters))


def genus(curve):
    """Return the smooth projective genus of a PlaneCurve."""
    return curve.genus


def normalization(curve, **limits):
    """Return exact plane-curve normalization analysis, including singularities and adjoints."""
    return curve.normalization(**limits)


def parametrize(curve, parameter=None, **options):
    """Return a rational affine parametrization of a genus-zero plane curve."""
    return curve.parametrize(parameter, **options)


__all__ = ["PlaneCurve", "genus", "normalization", "parametrize"]
