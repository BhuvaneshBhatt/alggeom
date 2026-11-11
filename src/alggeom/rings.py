"""Elements and coefficient-preserving homomorphisms of affine coordinate rings."""

from dataclasses import dataclass
from functools import cached_property

import sympy as sp
from algroots.domains import field_element
from sympy.polys.polyerrors import BasePolynomialError

from .ideals import CoordinateRing, Ideal, PolynomialRing


def _coordinate_ring(ring):
    if isinstance(ring, PolynomialRing):
        return ring.quotient()
    if not isinstance(ring, CoordinateRing):
        raise TypeError("a polynomial or coordinate ring is required")
    return ring


def _coefficient_embedding(source, target):
    """Check the specified embedding of rational, algebraic and transcendental constants."""
    if source.characteristic() != target.characteristic():
        raise ValueError("coefficient embeddings must preserve characteristic")
    ground = source.domain if source.is_FractionField else source
    from algroots.finite_fields import FiniteFieldExtension

    if isinstance(ground, FiniteFieldExtension) and ground != (
        target.domain if target.is_FractionField else target
    ):
        raise ValueError("finite-extension embeddings require the same declared field")
    generators = [ground.ext.as_expr()] if ground.is_AlgebraicField else []
    if getattr(source, "is_AlgebraicFunctionField", False):
        if not getattr(target, "is_AlgebraicFunctionField", False) or source != target:
            raise ValueError("algebraic function fields require an explicit coefficient embedding")
        generators.extend(source.parameters + source.tower_symbols)
    if source.is_FractionField:
        generators.extend(source.symbols)
    for generator in generators:
        field_element(generator, target)


@dataclass(frozen=True, eq=False)
class CoordinateElement:
    """An element of S/I with exact arithmetic and equality modulo I.

    Arithmetic requires the same parent coordinate ring, or an expression
    coercible into it. as_expr() removes the parent and returns a polynomial
    representative. Elements are unhashable because equal quotient rings may
    use different term orders and hence different representatives.
    """

    ring: CoordinateRing
    expression: sp.Expr

    __hash__ = None
    _op_priority = 1000

    def __post_init__(self):
        if not isinstance(self.ring, CoordinateRing):
            raise TypeError("an element requires a CoordinateRing parent")
        object.__setattr__(self, "expression", self.ring.normal_form(self.expression))

    def as_expr(self):
        """Return the reduced expression without its coordinate-ring parent."""
        return self.expression

    def __str__(self):
        return str(self.expression)

    def __bool__(self):
        return self.expression != 0

    def __eq__(self, other):
        from .rational_functions import RationalFunction

        if isinstance(other, RationalFunction):
            return other.parent.base == self.ring and other == self
        if isinstance(other, CoordinateElement):
            if self.ring != other.ring:
                return False
            other = other.expression
        try:
            return self.ring.normal_form(self.expression - other) == 0
        except (TypeError, ValueError, BasePolynomialError):
            return False

    def __add__(self, other):
        from .rational_functions import RationalFunction

        if isinstance(other, RationalFunction):
            return other + self
        other = self.ring.element(other)
        return self.ring.element(self.expression + other.expression)

    __radd__ = __add__

    def __neg__(self):
        return self.ring.element(-self.expression)

    def __sub__(self, other):
        from .rational_functions import RationalFunction

        if isinstance(other, RationalFunction):
            return self + (-other)
        return self + (-self.ring.element(other))

    def __rsub__(self, other):
        return self.ring.element(other) - self

    def __mul__(self, other):
        from .rational_functions import RationalFunction

        if isinstance(other, RationalFunction):
            return other * self
        other = self.ring.element(other)
        return self.ring.element(self.expression * other.expression)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, (int, sp.Integer)):
            raise TypeError("coordinate-ring powers require an integer exponent")
        if exponent < 0:
            return self.inverse() ** (-exponent)
        return self.ring.element(self.expression ** int(exponent))

    @cached_property
    def _inverse(self):
        ideal = self.ring.ideal + self.ring.ambient.ideal(self.expression)
        witness = ideal.membership(1)
        if witness is None:
            return None
        if self.ring.ideal.is_unit:
            return self.ring.zero
        # The appended element may coincide with a stored ideal generator.
        # That happens only for the zero class, which is not a unit here.
        return self.ring.element(witness[-1])

    @property
    def is_unit(self):
        """Whether an inverse exists in this coordinate ring, without assuming finiteness."""
        return self._inverse is not None

    def inverse(self):
        """Return the inverse in S/I, raising ZeroDivisionError for a nonunit."""
        if self._inverse is None:
            raise ZeroDivisionError("the coordinate-ring element is not a unit")
        return self._inverse

    def __truediv__(self, other):
        return self * self.ring.element(other).inverse()

    def __rtruediv__(self, other):
        return self.ring.element(other) * self.inverse()


@dataclass(frozen=True)
class RingHomomorphism:
    """A homomorphism S/I -> T/J defined by images of the S coordinates.

    Source and target have the same declared coefficient field. Construction
    checks all source relations. Ideal outputs are inverse-image ideals in the
    relevant ambient polynomial ring, including that ring's defining ideal.
    """

    source: CoordinateRing
    target: CoordinateRing
    images: tuple[CoordinateElement, ...]

    def __post_init__(self):
        object.__setattr__(self, "source", _coordinate_ring(self.source))
        object.__setattr__(self, "target", _coordinate_ring(self.target))
        if self.source.domain != self.target.domain:
            raise ValueError(
                "homomorphisms require the same coefficient field; use base_change first"
            )
        images = tuple(self.target.element(v) for v in self.images)
        if len(images) != len(self.source.variables):
            raise ValueError("one image is required for each source coordinate")
        object.__setattr__(self, "images", images)
        if any(self._substitute(g) != 0 for g in self.source.ideal.generators):
            raise ValueError("the images do not satisfy the source relations")

    @cached_property
    def _substitution(self):
        return dict(zip(self.source.variables, (v.as_expr() for v in self.images), strict=True))

    def _substitute(self, expression):
        # Inputs are polynomial representatives. Replacing coordinate symbols
        # simultaneously keeps overlapping source and target names independent.
        return self.target.element(sp.sympify(expression).xreplace(self._substitution))

    def __call__(self, element):
        return self._substitute(self.source.element(element).as_expr())

    @cached_property
    def _source_symbols(self):
        return tuple(sp.Dummy("source_" + v.name) for v in self.source.variables)

    @cached_property
    def graph_ideal(self):
        """The graph in target-first lex order, with fresh source coordinate symbols."""
        ring = PolynomialRing(
            (*self.target.variables, *self._source_symbols), self.target.domain, "lex"
        )
        renamed = dict(zip(self.source.variables, self._source_symbols, strict=True))
        return ring.ideal(
            *self.target.ideal.generators,
            *(g.xreplace(renamed) for g in self.source.ideal.generators),
            *(
                v - image.as_expr()
                for v, image in zip(self._source_symbols, self.images, strict=True)
            ),
        )

    def _contract(self, ideal):
        result = ideal.eliminate(*self.target.variables)
        renamed = dict(zip(self._source_symbols, self.source.variables, strict=True))
        return self.source.ambient.ideal(*(g.xreplace(renamed) for g in result.generators))

    @cached_property
    def kernel(self):
        """The kernel as an ambient ideal containing the source defining ideal."""
        return self._contract(self.graph_ideal)

    def preimage(self, ideal):
        """Contract a target ideal along this homomorphism."""
        if (
            not isinstance(ideal, Ideal)
            or ideal.ring.variables != self.target.variables
            or ideal.ring.domain != self.target.domain
        ):
            raise ValueError("the ideal must belong to the target ambient ring")
        return self._contract(self.graph_ideal + ideal.change_ring(self.graph_ideal.ring))

    def extend(self, ideal):
        """Generate the target ideal from images of a source ideal's generators."""
        if (
            not isinstance(ideal, Ideal)
            or ideal.ring.variables != self.source.variables
            or ideal.ring.domain != self.source.domain
        ):
            raise ValueError("the ideal must belong to the source ambient ring")
        return self.target.ideal + self.target.ambient.ideal(
            *(self(g).as_expr() for g in ideal.generators)
        )

    @property
    def is_injective(self):
        """Whether the map is injective, using exact kernel contraction."""
        return self.kernel == self.source.ideal

    def lift(self, element):
        """Return one polynomial preimage of a target element, or raise ValueError."""
        reduced = self.graph_ideal.normal_form(self.target.element(element).as_expr())
        if reduced.free_symbols & set(self.target.variables):
            raise ValueError("the element does not belong to the homomorphism image")
        return self.source.element(
            reduced.xreplace(dict(zip(self._source_symbols, self.source.variables, strict=True)))
        )

    @property
    def is_surjective(self):
        """Whether all target coordinate classes have polynomial preimages."""
        try:
            for generator in self.target.gens:
                self.lift(generator)
            return True
        except ValueError:
            return False

    def inverse(self):
        """Return the inverse of a bijective homomorphism."""
        if not self.is_injective:
            raise ValueError("the homomorphism is not injective")
        return RingHomomorphism(
            self.target, self.source, tuple(self.lift(g) for g in self.target.gens)
        )

    def compose(self, inner):
        """Return self after inner, requiring the intermediate coordinate rings to agree."""
        if not isinstance(inner, RingHomomorphism) or inner.target != self.source:
            raise ValueError("the intermediate coordinate rings must agree")
        return RingHomomorphism(inner.source, self.target, tuple(self(v) for v in inner.images))

    def base_change(self, domain):
        """Extend source, target and coordinate images to the same exact coefficient field."""
        return RingHomomorphism(
            self.source.base_change(domain),
            self.target.base_change(domain),
            tuple(v.as_expr() for v in self.images),
        )


__all__ = ["CoordinateElement", "RingHomomorphism"]
