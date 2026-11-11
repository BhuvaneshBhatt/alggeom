"""Exact arithmetic on principal opens of arbitrary parameter quotient rings.

For A = QQ[a]/I and h = product(inverted), the localization A[h^-1]
has kernel (I:h^infinity)/I. Consequently equality is decidable by reducing
cross products modulo the *saturated* ideal. No primality, reducedness, or
non-zero-divisor hypothesis is imposed: components on which h vanishes
identically disappear upon localization, as algebraic geometry requires.

This is a general, finitely-generated localization, NOT the total fraction
ring of a non-domain, and it does not perform a comprehensive cover search.
"""

from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .ideals import Ideal


@dataclass(frozen=True)
class ParameterLocalization:
    """QQ[a]/I localized at finitely many explicitly inverted polynomials."""

    ideal: Ideal
    inverted: tuple[sp.Expr, ...] = ()

    def __post_init__(self):
        if not isinstance(self.ideal, Ideal) or self.ideal.ring.domain != sp.QQ:
            raise TypeError("localization requires a polynomial ideal over QQ")
        if not self.ideal.ring.variables:
            raise ValueError("at least one parameter coordinate is required")
        values = tuple(self.ideal.ring.convert(value) for value in self.inverted)
        if any(value == 0 for value in values):
            raise ZeroDivisionError("cannot invert the zero polynomial")
        object.__setattr__(self, "inverted", values)
        if self.kernel_ideal.is_unit:
            raise ZeroDivisionError("the specified principal open is empty")

    @cached_property
    def kernel_ideal(self):
        """Preimage of zero under A -> S^-1 A, computed by exact saturation."""
        if not self.inverted:
            return self.ideal
        product = sp.prod(self.inverted)
        return self.ideal.saturate(product)

    def reduce(self, polynomial):
        """Compute the canonical class modulo the localization kernel."""
        return sp.expand(self.kernel_ideal.normal_form(self.ideal.ring.convert(polynomial)))

    def element(self, numerator, powers=None):
        """Represent numerator / product(inverted[i]**powers[i])."""
        if powers is None:
            powers = (0,) * len(self.inverted)
        powers = tuple(powers)
        if len(powers) != len(self.inverted) or any(
            not isinstance(power, int) or power < 0 for power in powers
        ):
            raise ValueError("one nonnegative integer power per inverted generator is required")
        return LocalizedParameterElement(self, self.reduce(numerator), powers)

    def coerce(self, value):
        """Coerce a polynomial or an element of this exact localization."""
        if isinstance(value, LocalizedParameterElement):
            if value.localization != self:
                raise TypeError("different localizations require an explicit restriction map")
            return value
        return self.element(value)

    def inverse_generator(self, index):
        """Return the certified inverse of one declared localization generator."""
        if not isinstance(index, int) or not 0 <= index < len(self.inverted):
            raise IndexError("unknown localization generator")
        return self.element(1, tuple(int(j == index) for j in range(len(self.inverted))))

    def localize_at(self, *polynomials):
        """Refine to the intersection with additional principal opens."""
        return ParameterLocalization(self.ideal, (*self.inverted, *polynomials))

    @property
    def condition(self):
        """Constructible condition of this quotient's principal open."""
        return sp.And(
            *(sp.Eq(value, 0) for value in self.ideal.groebner_basis),
            *(sp.Ne(value, 0) for value in self.inverted),
        )

    def split_at(self, polynomial):
        """Partition this locus into D(polynomial) and V(polynomial).

        Each branch is a genuine localization or quotient, even when the
        parameter ideal has embedded components or zero divisors. An empty
        branch is represented by None; coverage/disjointness follows from
        f != 0 or f == 0 on every geometric parameter point.
        """
        polynomial = self.ideal.ring.convert(polynomial)
        try:
            opened = self.localize_at(polynomial)
        except ZeroDivisionError:
            opened = None
        try:
            closed = ParameterLocalization(
                self.ideal + self.ideal.ring.ideal(polynomial), self.inverted
            )
        except ZeroDivisionError:
            closed = None
        return opened, closed

    def split_ideal(self, exceptional):
        """Disjointly split D(J) and V(J) for any finitely generated ideal J.

        D(J) = D(f_1) union [V(f_1) intersect D(f_2)] union ...;
        the pieces are disjoint by construction. The last branch is V(J),
        without assuming J is radical, prime or equidimensional. This
        provides the exact branching infrastructure, not a Gröbner cover on
        zero-divisor quotient rings by itself.
        """
        if not isinstance(exceptional, Ideal):
            exceptional = self.ideal.ring.ideal(exceptional)
        if (
            exceptional.ring.variables != self.ideal.ring.variables
            or exceptional.ring.domain != self.ideal.ring.domain
        ):
            raise ValueError("exceptional ideal belongs to another parameter ring")
        opened_branches = []
        closed = self
        for generator in exceptional.generators:
            if closed is None:
                break
            opened, closed = closed.split_at(generator)
            if opened is not None:
                opened_branches.append(opened)
        return tuple(opened_branches), closed

    def restrict(self, element):
        """Transport a fraction from an earlier prefix localization to this one."""
        if not isinstance(element, LocalizedParameterElement):
            return self.element(element)
        old = element.localization
        if old.ideal != self.ideal or self.inverted[: len(old.inverted)] != old.inverted:
            raise TypeError("restriction requires the same quotient and a prefix open")
        extra = (0,) * (len(self.inverted) - len(old.inverted))
        return self.element(element.numerator, (*element.powers, *extra))


@dataclass(frozen=True, eq=False)
class LocalizedParameterElement:
    """A fraction represented by a numerator and explicit inversion powers."""

    localization: ParameterLocalization
    numerator: sp.Expr
    powers: tuple[int, ...]

    def _denominator(self):
        return sp.prod(
            value**power
            for value, power in zip(self.localization.inverted, self.powers, strict=True)
        )

    def _scaled_numerator(self, target_powers):
        return self.numerator * sp.prod(
            value ** (target - existing)
            for value, existing, target in zip(
                self.localization.inverted, self.powers, target_powers, strict=True
            )
        )

    def __add__(self, other):
        other = self.localization.coerce(other)
        powers = tuple(max(a, b) for a, b in zip(self.powers, other.powers, strict=True))
        return self.localization.element(
            self._scaled_numerator(powers) + other._scaled_numerator(powers), powers
        )

    __radd__ = __add__

    def __neg__(self):
        return self.localization.element(-self.numerator, self.powers)

    def __sub__(self, other):
        return self + -self.localization.coerce(other)

    def __rsub__(self, other):
        return self.localization.coerce(other) - self

    def __mul__(self, other):
        other = self.localization.coerce(other)
        powers = tuple(a + b for a, b in zip(self.powers, other.powers, strict=True))
        return self.localization.element(self.numerator * other.numerator, powers)

    __rmul__ = __mul__

    def __eq__(self, other):
        try:
            other = self.localization.coerce(other)
        except (TypeError, ValueError):
            return NotImplemented
        return (
            self.localization.reduce(
                self.numerator * other._denominator() - other.numerator * self._denominator()
            )
            == 0
        )

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("only nonnegative integer powers are supported; invert explicitly")
        return self.localization.element(
            self.numerator**power, tuple(p * power for p in self.powers)
        )

    def is_zero(self):
        """Decide whether the fraction is zero on this localized quotient."""
        return self.localization.reduce(self.numerator) == 0

    def as_expr(self):
        """Return a rational-expression representative; equality uses saturation."""
        return sp.cancel(self.numerator / self._denominator())


__all__ = ["LocalizedParameterElement", "ParameterLocalization"]
