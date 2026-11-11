"""Parent-aware rational functions, principal localizations and local rings."""

from dataclasses import dataclass
from functools import cached_property

import sympy as sp
from sympy.polys.polyerrors import BasePolynomialError

from ._algebra import integer
from .ideals import CoordinateRing, Ideal, PolynomialRing
from .rings import CoordinateElement


def _base(value):
    if isinstance(value, PolynomialRing):
        return value.quotient()
    if isinstance(value, Ideal):
        return CoordinateRing(value)
    if not isinstance(value, CoordinateRing):
        raise TypeError("a polynomial or coordinate ring is required")
    return value


class _FractionParent:
    def _expression(self, value):
        if isinstance(value, CoordinateElement):
            if value.ring is not self.base and value.ring != self.base:
                raise TypeError("coordinate elements require the same base coordinate ring")
            return value.as_expr()
        if isinstance(value, RationalFunction) and value.parent.base != self.base:
            raise TypeError("rational functions require the same base coordinate ring")
        return value.as_expr() if hasattr(value, "as_expr") else value

    def fraction(self, numerator, denominator=1):
        """Construct a fraction, checking its base ring and allowed denominators.

        Coordinate elements must belong to the base coordinate ring. Use
        as_expr() explicitly to reinterpret a representative in another ring.
        """
        denominator = self._expression(denominator)
        if isinstance(numerator, RationalFunction):
            if numerator.parent.base != self.base:
                raise TypeError("rational functions require the same base coordinate ring")
            denominator = sp.sympify(denominator) * numerator.denominator
            numerator = numerator.numerator
        numerator = self._expression(numerator)
        return RationalFunction(self, numerator, denominator)

    def __call__(self, expression):
        if isinstance(expression, RationalFunction):
            if expression.parent is self:
                return expression
            return self.fraction(expression)
        if isinstance(expression, CoordinateElement):
            return self.fraction(expression)
        expression = (
            expression.as_expr() if hasattr(expression, "as_expr") else sp.sympify(expression)
        )
        numerator, denominator = sp.together(expression).as_numer_denom()
        return self.fraction(numerator, denominator)

    @property
    def gens(self):
        """The images of the base coordinate generators."""
        return tuple(self(v) for v in self.base.variables)

    @property
    def zero(self):
        """The additive identity."""
        return self(0)

    @property
    def one(self):
        """The multiplicative identity."""
        return self(1)

    @property
    def domain(self):
        """The coefficient field."""
        return self.base.domain

    def _normalize(self, numerator, denominator):
        numerator = self.kernel.normal_form(numerator)
        denominator = self.kernel.normal_form(denominator)
        if self.kernel.is_unit:
            return sp.S.Zero, sp.S.One
        if denominator == 0:
            raise ZeroDivisionError("the denominator is zero in the localized ring")
        return numerator, denominator

    def _cancel(self, numerator, denominator):
        ambient = self.base.ambient
        left, right = ambient._poly(numerator), ambient._poly(denominator)
        common = left.gcd(right)
        left, right = left.exquo(common), right.exquo(common)
        inverse = ambient._poly_ring.domain.exquo(ambient._poly_ring.domain.one, right.LC)
        return left.mul_ground(inverse).as_expr(), right.mul_ground(inverse).as_expr()

    def _regular_pair(self, numerator, denominator):
        """Find a regular presentation through the denominator ideal (b):(a).

        A rational function can be regular even when its displayed denominator
        vanishes at the point. A witness u*a=b*v with u invertible in the
        localization supplies the equivalent fraction v/u. Zero-divisor
        denominators are rejected because that identity would not determine
        a unique rational function.
        """
        if self._is_inverted(denominator):
            return numerator, denominator
        if self.kernel.quotient(denominator) != self.kernel:
            raise ZeroDivisionError(
                "a zero-divisor denominator does not define a rational function"
            )
        multiples = self.kernel + self.base.ambient.ideal(denominator)
        denominator_ideal = multiples.quotient(numerator)
        witness = self._unit_in(denominator_ideal)
        if witness is None:
            raise ValueError("the rational function is not regular in this localization")
        coefficients = multiples.membership(witness * numerator)
        return coefficients[-1], witness


@dataclass(frozen=True)
class FractionField(_FractionParent):
    """The field of rational functions of a nonzero integral coordinate ring.

    Fractions are compared by cross multiplication modulo the defining ideal.
    Polynomial gcd cancellation is only a representation simplification; a
    unique numerator/denominator form is not assumed in a nonfactorial ring.
    """

    base: CoordinateRing

    def __post_init__(self):
        object.__setattr__(self, "base", _base(self.base))
        if not self.base.is_integral:
            raise ValueError("a fraction field requires an integral coordinate ring")

    @property
    def kernel(self):
        """The zero ideal of the base, lifted to its polynomial ambient ring."""
        return self.base.ideal

    def _regular_pair(self, numerator, denominator):
        return numerator, denominator


@dataclass(frozen=True)
class Localization(_FractionParent):
    """The ring obtained by inverting finitely many functions in a coordinate ring.

    An explicit affine model adjoins u with u*product(denominators)=1. Kernel
    saturation preserves correct equality when localization removes components.
    """

    base: CoordinateRing
    denominators: tuple[sp.Expr, ...]
    max_power: int = 256

    def __post_init__(self):
        object.__setattr__(self, "base", _base(self.base))
        object.__setattr__(
            self, "denominators", tuple(self.base.element(v).as_expr() for v in self.denominators)
        )
        object.__setattr__(self, "max_power", integer(self.max_power, "max_power", 1))

    @property
    def inverted(self):
        """The product of the specified invertible functions."""
        return sp.prod(self.denominators)

    @cached_property
    def kernel(self):
        """The contraction of zero from the localization: I:f^infinity."""
        return self.base.ideal.saturate(self.inverted)

    @cached_property
    def affine_model(self):
        """An affine coordinate-ring presentation of the localization."""
        u = sp.Dummy("inverse")
        ring = PolynomialRing((*self.base.variables, u), self.base.domain)
        return ring.quotient(*self.base.ideal.generators, u * self.inverted - 1)

    @property
    def natural_map(self):
        """The canonical map from the base ring to the affine localization model."""
        return self.base.hom(self.base.variables, self.affine_model)

    def _is_inverted(self, denominator):
        return self.affine_model(denominator).is_unit

    def _unit_in(self, ideal):
        if self.affine_model.ideal.is_unit:
            return sp.S.One
        extended = self.affine_model.ideal + self.affine_model.ambient.ideal(*ideal.generators)
        if not extended.is_unit:
            return None
        power = 0
        while power <= self.max_power:
            candidate = self.inverted**power
            if ideal.contains(candidate):
                return candidate
            power = (
                1
                if power == 0
                else min(2 * power, self.max_power)
                if power < self.max_power
                else self.max_power + 1
            )
        raise RuntimeError("localization witness exceeded max_power")


@dataclass(frozen=True)
class LocalRing(_FractionParent):
    """The localization of a coordinate ring at a specified prime ideal.

    The prime is represented in the ambient polynomial ring and must contain
    the defining ideal. Component contraction computes equality with nilpotents
    and zero divisors. General local standard bases and completion are not
    provided by this arithmetic interface.
    """

    base: CoordinateRing
    prime: Ideal

    def __post_init__(self):
        from .components import is_prime

        object.__setattr__(self, "base", _base(self.base))
        if not isinstance(self.prime, Ideal) or self.prime.ring != self.base.ambient:
            raise ValueError("a prime in the base ambient ring is required")
        if not self.base.ideal.is_subset(self.prime) or not is_prime(self.prime):
            raise ValueError("the prime must contain the base defining ideal")

    @cached_property
    def kernel(self):
        """The contraction of zero after localizing at the prime."""
        from .components import primary_decomposition

        result = self.base.ambient.ideal(1)
        for component in primary_decomposition(self.base.ideal):
            if component.radical().is_subset(self.prime):
                result = result.intersection(component)
        return result

    @property
    def maximal_ideal(self):
        """The ambient inverse image of the maximal ideal of this local ring."""
        return self.prime

    @property
    def residue_field(self):
        """The fraction field of the integral residue coordinate ring."""
        return FractionField(CoordinateRing(self.prime))

    def _is_inverted(self, denominator):
        return not self.prime.contains(denominator)

    def _unit_in(self, ideal):
        return next((g for g in ideal.generators if not self.prime.contains(g)), None)


@dataclass(frozen=True, eq=False)
class RationalFunction:
    """A fraction in an explicit fraction field, localization or local ring."""

    parent: _FractionParent
    numerator: sp.Expr
    denominator: sp.Expr = sp.S.One

    __hash__ = None
    _op_priority = 1000

    def __post_init__(self):
        if not isinstance(self.parent, _FractionParent):
            raise TypeError("a fraction-field or localization parent is required")
        numerator, denominator = self.parent._normalize(self.numerator, self.denominator)
        numerator, denominator = self.parent._regular_pair(numerator, denominator)
        numerator, denominator = self.parent._normalize(numerator, denominator)
        numerator, denominator = self.parent._cancel(numerator, denominator)
        object.__setattr__(self, "numerator", numerator)
        object.__setattr__(self, "denominator", denominator)

    def as_expr(self):
        """Return an ambient rational expression, removing the quotient-ring parent."""
        return self.numerator / self.denominator

    def __bool__(self):
        return not self.parent.kernel.contains(self.numerator)

    def __eq__(self, other):
        if isinstance(other, RationalFunction) and other.parent != self.parent:
            return False
        try:
            other = self.parent(other)
        except (TypeError, ValueError, ZeroDivisionError, BasePolynomialError):
            return False
        return self.parent.kernel.contains(
            self.numerator * other.denominator - other.numerator * self.denominator
        )

    def __add__(self, other):
        other = self.parent(other)
        return self.parent.fraction(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )

    __radd__ = __add__

    def __neg__(self):
        return self.parent.fraction(-self.numerator, self.denominator)

    def __sub__(self, other):
        return self + (-self.parent(other))

    def __rsub__(self, other):
        return self.parent(other) - self

    def __mul__(self, other):
        other = self.parent(other)
        return self.parent.fraction(
            self.numerator * other.numerator, self.denominator * other.denominator
        )

    __rmul__ = __mul__

    def inverse(self):
        """Invert a unit; a nonunit or zero raises ZeroDivisionError."""
        try:
            return self.parent.fraction(self.denominator, self.numerator)
        except ValueError as exc:
            raise ZeroDivisionError("the element is not a unit in this localization") from exc

    @property
    def is_unit(self):
        """Whether the inverse belongs to the same parent."""
        try:
            self.inverse()
        except ZeroDivisionError:
            return False
        return True

    def __truediv__(self, other):
        return self * self.parent(other).inverse()

    def __rtruediv__(self, other):
        return self.parent(other) * self.inverse()

    def __pow__(self, exponent):
        exponent = integer(exponent, "exponent", minimum=None)
        if exponent < 0:
            return self.inverse() ** (-exponent)
        return self.parent.fraction(self.numerator**exponent, self.denominator**exponent)

    def to_coordinate(self):
        """Return a base-ring element if this rational function is globally regular."""
        ambient = self.parent.base.ambient
        multiples = self.parent.base.ideal + ambient.ideal(self.denominator)
        witness = multiples.membership(self.numerator)
        if witness is None:
            raise ValueError("the rational function is not represented in the base coordinate ring")
        if self.parent.base.ideal.contains(self.denominator):
            raise ValueError("the displayed denominator vanishes in the base ring")
        return self.parent.base(witness[-1])


def fraction_field(value):
    """Return the rational-function field of an integral coordinate ring."""
    return FractionField(_base(value))


def local_ring(value, point_or_prime):
    """Return a local ring at a prime or at an exact affine/projective point."""
    if isinstance(point_or_prime, Ideal):
        return LocalRing(_base(value), point_or_prime)
    from .local import _affine_point

    ideal, values = _affine_point(value, point_or_prime)
    prime = ideal.ring.ideal(*(v - p for v, p in zip(ideal.ring.variables, values, strict=True)))
    return LocalRing(CoordinateRing(ideal), prime)


__all__ = [
    "FractionField",
    "Localization",
    "LocalRing",
    "RationalFunction",
    "fraction_field",
    "local_ring",
]
