"""Exact polynomial ideals and coordinate rings over supported exact coefficient fields."""

from dataclasses import dataclass
from functools import cached_property

import sympy as sp
from algroots.domains import field_element, polynomial, polynomial_domain
from algroots.quotient import QuotientAlgebra
from sympy.polys.groebnertools import groebner
from sympy.polys.rings import PolyRing

from ._algebra import coordinates, exact_field


@dataclass(frozen=True)
class PolynomialRing:
    """An ordered polynomial ring over a supported exact coefficient field.

    Variables may be empty, representing the coefficient field itself. The term
    order controls normal forms; it is not part of mathematical ideal equality.
    Coefficient parameters and polynomial variables must be disjoint.
    """

    variables: tuple[sp.Symbol, ...]
    domain: object = sp.QQ
    order: str = "grevlex"

    def __post_init__(self):
        object.__setattr__(self, "variables", coordinates(self.variables))
        object.__setattr__(self, "domain", exact_field(self.domain))
        if self.order not in {"lex", "grlex", "grevlex"}:
            raise ValueError("order must be lex, grlex, or grevlex")
        parameters = self.domain.symbols if self.domain.is_FractionField else ()
        parameters = (
            *parameters,
            *getattr(self.domain, "parameters", ()),
            *getattr(self.domain, "tower_symbols", ()),
        )
        ground = self.domain.domain if self.domain.is_FractionField else self.domain
        if getattr(ground, "is_FiniteExtension", False):
            parameters = (*parameters, ground.symbol)
        if any(v == p or v.name == p.name for v in self.variables for p in parameters):
            raise ValueError("coefficient parameters and polynomial coordinates must be distinct")

    @cached_property
    def _poly_ring(self):
        return PolyRing(self.variables, polynomial_domain(self.domain), self.order)

    def _poly(self, expression):
        if not self.variables:
            return self._poly_ring.ground_new(field_element(expression, self._poly_ring.domain))
        expression = (
            expression.as_expr() if isinstance(expression, sp.Poly) else sp.sympify(expression)
        )
        if expression.has(sp.Float):
            raise ValueError("exact polynomial coefficients are required")
        return self._poly_ring.from_dict(
            polynomial(expression, self.variables, self._poly_ring.domain).rep.to_dict()
        )

    def convert(self, expression):
        """Return the canonical exact expression, rejecting poles in the ring variables."""
        return self._poly(expression).as_expr()

    def ideal(self, *generators, groebner_backend="native", groebner_limits=None):
        """Generate an ideal with an explicit verified Gröbner construction policy.

        Accept expressions or one list/tuple. Backend selection is stored on
        this ideal and controls its cached basis and normal forms.
        """
        if len(generators) == 1 and isinstance(generators[0], (tuple, list)):
            generators = tuple(generators[0])
        return Ideal(self, tuple(generators), groebner_backend, groebner_limits)

    def quotient(self, *generators, groebner_backend="native", groebner_limits=None):
        """Return a coordinate ring with the selected verified basis backend."""
        return CoordinateRing(
            self.ideal(
                *generators, groebner_backend=groebner_backend, groebner_limits=groebner_limits
            )
        )

    def with_order(self, order):
        """Use another term order with the same coefficient field and coordinates."""
        return PolynomialRing(self.variables, self.domain, order)


@dataclass(frozen=True, eq=False)
class Ideal:
    """A finitely generated ideal, with mathematical equality and exact normal forms.

    Zero generators are omitted; the empty generator tuple is the zero ideal.
    Membership witnesses refer to the stored generators, not only a computed
    Groebner basis. Equality requires identical coordinates and coefficients,
    but permits different term orders.
    """

    ring: PolynomialRing
    generators: tuple[sp.Expr, ...] = ()
    groebner_backend: str = "native"
    groebner_limits: object = None

    __hash__ = None

    def __post_init__(self):
        from algroots.basis import _validate_backend

        _validate_backend(self.groebner_backend, self.groebner_limits)
        if not isinstance(self.ring, PolynomialRing):
            raise TypeError("an ideal requires a PolynomialRing")
        values = tuple(self.ring.convert(g) for g in self.generators)
        object.__setattr__(self, "generators", tuple(dict.fromkeys(g for g in values if g != 0)))

    @cached_property
    def _basis(self):
        if not self.ring.variables:
            return (self.ring._poly_ring.one,) if self.generators else ()
        if self.groebner_backend != "sympy" or self.groebner_limits is not None:
            return tuple(
                self.ring._poly(p)
                for p in self.groebner_basis_from(
                    self.groebner_backend, limits=self.groebner_limits
                )
            )
        inputs = [self.ring._poly(g) for g in self.generators]
        return tuple(groebner(inputs, self.ring._poly_ring))

    def groebner_basis_from(self, backend="native", *, limits=None):
        """Return a reduced basis constructed and verified through algroots.

        Native requests over a prime field run directly in that characteristic.
        Rational requests use modular lifting. Every accelerated result passes
        an independent exact check before it enters geometric calculations.
        """
        from algroots import groebner_basis

        basis = groebner_basis(
            self.generators,
            self.ring.variables,
            domain=self.ring.domain,
            order=str(self.ring._poly_ring.order),
            backend=backend,
            limits=limits,
        )
        return tuple(poly.as_expr() for poly in basis.polys)

    def groebner_basis_proof_from(
        self,
        backend="native",
        *,
        limits=None,
        max_witness_basis=128,
        max_witness_pairs=500,
    ):
        """Return a basis with independently checked generator-combination witnesses.

        Algroots constructs the candidate and its membership witnesses.
        Exact replay in alggeom establishes both ideal inclusions before the proof is returned. Use
        ``backend="sympy"`` for a portable candidate calculation.
        """
        from algroots.groebner import compute_groebner_basis

        from .groebner_certificates import (
            GroebnerBasisProof,
            construct_membership_witnesses,
        )

        if not self.ring.variables:
            raise ValueError("the backend requires polynomial variables")
        original = tuple(self.ring._poly(g) for g in self.generators)
        candidate = compute_groebner_basis(
            original, self.ring._poly_ring, backend=backend, limits=limits, certify=False
        )
        witnesses = construct_membership_witnesses(
            original,
            candidate,
            self.ring._poly_ring,
            max_basis=max_witness_basis,
            max_pairs=max_witness_pairs,
        )
        proof = GroebnerBasisProof(
            tuple(g.as_expr() for g in candidate),
            tuple(
                tuple(
                    coefficient.as_expr() if hasattr(coefficient, "as_expr") else dict(coefficient)
                    for coefficient in row
                )
                for row in witnesses
            ),
        )
        if not proof.verify(original, self.ring._poly_ring):
            raise ArithmeticError(f"{backend} basis failed independent proof verification")
        return proof

    @property
    def groebner_basis(self):
        """The reduced monic Groebner basis as a tuple of expressions."""
        return tuple(p.as_expr() for p in self._basis)

    @property
    def is_zero(self):
        """Whether the ideal is zero."""
        return not self._basis

    @property
    def is_unit(self):
        """Whether the ideal is the entire polynomial ring."""
        return len(self._basis) == 1 and self._basis[0] == self.ring._poly_ring.one

    @cached_property
    def krull_dimension(self):
        """The dimension of the quotient S/I; -1 denotes the zero ring."""
        from ._monomials import dimension

        return dimension((p.LM for p in self._basis), len(self.ring.variables))

    def divide(self, expression):
        """Return quotients against groebner_basis and the unique reduced remainder."""
        quotients, remainder = self.ring._poly(expression).div(list(self._basis))
        return tuple(q.as_expr() for q in quotients), remainder.as_expr()

    def normal_form(self, expression):
        """Return the reduced polynomial representative modulo this ideal."""
        return self.ring._poly(expression).rem(self._basis).as_expr()

    def contains(self, expression):
        """Decide exact ideal membership; invalid coefficient/ring inputs raise."""
        return self.normal_form(expression) == 0

    def is_subset(self, other):
        """Decide containment in another ideal of the same ambient ring."""
        other = self._coerce(other)
        return all(other.contains(g) for g in self.generators)

    def radical_contains(self, expression):
        """Decide f in radical(I) by the exact saturation criterion I:f^infinity=S."""
        return self.saturate(expression).is_unit

    def radical(self, **limits):
        """Return the radical over the declared coefficient field."""
        from .components import radical

        return radical(self, **limits)

    def minimal_primes(self, **limits):
        """Return the prime ideals of the reduced irreducible components."""
        from .components import minimal_primes

        return minimal_primes(self, **limits)

    def associated_primes(self, **limits):
        """Return all associated primes, including embedded primes."""
        from .components import associated_primes

        return associated_primes(self, **limits)

    def primary_decomposition(self, **limits):
        """Return an irredundant primary decomposition retaining embedded components."""
        from .components import primary_decomposition

        return primary_decomposition(self, **limits)

    @cached_property
    def _module(self):
        base = self.ring._poly_ring.domain.old_poly_ring(
            *self.ring.variables, order=self.ring.order
        )
        return base.free_module(1).submodule(
            *[
                [base.convert(polynomial(g, self.ring.variables, base.dom).rep)]
                for g in self.generators
            ]
        )

    def membership(self, expression):
        """Return coefficients expressing membership in the stored generators, or None."""
        expression = self.ring.convert(expression)
        if not self.contains(expression):
            return None
        if not self.generators:
            return ()
        if not self.ring.variables:
            coefficient = field_element(expression, self.ring.domain) / field_element(
                self.generators[0], self.ring.domain
            )
            return (
                self.ring.domain.to_sympy(coefficient),
                *(sp.S.Zero for _ in self.generators[1:]),
            )
        base = self._module.ring
        value = base.convert(polynomial(expression, self.ring.variables, base.dom).rep)
        coefficients = self._module.in_terms_of_generators([value])
        return tuple(base.to_sympy(c) for c in coefficients)

    def _coerce(self, other):
        if not isinstance(other, Ideal):
            other = self._generated_ideal((other,))
        if self.ring.variables != other.ring.variables or self.ring.domain != other.ring.domain:
            raise ValueError("ideal operations require the same coefficient field and coordinates")
        return other

    def _generated_ideal(self, generators, ring=None):
        """Preserve construction policy through derived ideal operations."""
        return Ideal(
            self.ring if ring is None else ring,
            tuple(generators),
            self.groebner_backend,
            self.groebner_limits,
        )

    def change_ring(self, ring):
        """Embed the generators in another explicitly supplied polynomial ring."""
        return self._generated_ideal(self.generators, ring)

    def __eq__(self, other):
        if not isinstance(other, Ideal):
            return NotImplemented
        if self.ring.variables != other.ring.variables or self.ring.domain != other.ring.domain:
            return False
        return all(self.contains(g) for g in other.generators) and all(
            other.contains(g) for g in self.generators
        )

    def __contains__(self, expression):
        return self.contains(expression)

    def __add__(self, other):
        other = self._coerce(other)
        return self._generated_ideal((*self.generators, *other.generators))

    def __mul__(self, other):
        other = self._coerce(other)
        return self._generated_ideal(f * g for f in self.generators for g in other.generators)

    def eliminate(self, *variables):
        """Contract to the ring of retained coordinates by lexicographic elimination.

        Eliminating every coordinate returns an ideal in the coefficient
        field. Elimination describes a Zariski closure, not an exact image of
        points with inequalities or pole exclusions.
        """
        if len(variables) == 1 and isinstance(variables[0], (tuple, list)):
            variables = tuple(variables[0])
        variables = coordinates(variables)
        if not set(variables) <= set(self.ring.variables):
            raise ValueError("eliminated coordinates must belong to the ambient ring")
        keep = tuple(v for v in self.ring.variables if v not in variables)
        ordered = PolynomialRing((*variables, *keep), self.ring.domain, "lex")
        basis = self.change_ring(ordered)._basis
        expressions = [
            g.as_expr() for g in basis if all(not any(m[: len(variables)]) for m in g.monoms())
        ]
        return self._generated_ideal(
            expressions, PolynomialRing(keep, self.ring.domain, self.ring.order)
        )

    def intersection(self, other):
        """Return the exact intersection using an auxiliary-coordinate elimination."""
        other = self._coerce(other)
        t = sp.Dummy("intersection_coordinate")
        ring = PolynomialRing((t, *self.ring.variables), self.ring.domain, "lex")
        combined = self._generated_ideal(
            (*(t * g for g in self.generators), *((1 - t) * g for g in other.generators)), ring
        )
        return combined.eliminate(t).change_ring(self.ring)

    def quotient(self, other):
        """Return I:J; a zero divisor ideal J=0 gives the unit ideal."""
        other = self._coerce(other)
        result = self._generated_ideal((1,))
        for generator in other.generators:
            common = self.intersection(self._generated_ideal((generator,)))
            divisor = self.ring._poly(generator)
            principal = self._generated_ideal(g.exquo(divisor).as_expr() for g in common._basis)
            result = result.intersection(principal)
        return result

    def saturate(self, other):
        """Return I:J^infinity, supporting a polynomial or a finitely generated ideal.

        Principal saturation eliminates t from I+(1-t*f). For several
        generators intersect the principal saturations; saturating by their
        product would compute a different operation.
        """
        other = self._coerce(other)
        result = None
        for generator in other.generators:
            t = sp.Dummy("saturation_coordinate")
            ring = PolynomialRing((t, *self.ring.variables), self.ring.domain, "lex")
            principal = (
                self._generated_ideal((*self.generators, 1 - t * generator), ring)
                .eliminate(t)
                .change_ring(self.ring)
            )
            result = principal if result is None else result.intersection(principal)
        return self._generated_ideal((1,)) if result is None else result


@dataclass(frozen=True)
class CoordinateRing:
    """The possibly positive-dimensional algebra S/I.

    krull_dimension measures chains of prime ideals. vector_space_dimension
    is the length over the coefficient field and is infinite for a nonempty
    positive-dimensional algebra. The zero ring has Krull dimension -1 here.
    """

    ideal: Ideal

    def __post_init__(self):
        if not isinstance(self.ideal, Ideal):
            raise TypeError("a coordinate ring requires an Ideal")

    @property
    def ambient(self):
        """The ambient polynomial ring S."""
        return self.ideal.ring

    @property
    def variables(self):
        """The ordered ambient coordinate symbols."""
        return self.ambient.variables

    @property
    def domain(self):
        """The coefficient field, independent of the quotient's Krull dimension."""
        return self.ambient.domain

    def normal_form(self, expression):
        """Return the canonical representative of a coordinate-ring element."""
        return self.ideal.normal_form(expression)

    def element(self, expression):
        """Construct an element with this parent and a reduced representative."""
        from .rings import CoordinateElement

        if isinstance(expression, CoordinateElement):
            if expression.ring is self:
                return expression
            if expression.ring != self:
                raise TypeError("coordinate elements must have the same parent ring")
            expression = expression.as_expr()
        return CoordinateElement(self, expression)

    def __call__(self, expression):
        return self.element(expression)

    @cached_property
    def gens(self):
        """The residue classes of the ordered ambient coordinates."""
        return tuple(self.element(v) for v in self.variables)

    @property
    def zero(self):
        """The additive identity, including in the zero ring."""
        return self.element(0)

    @property
    def one(self):
        """The multiplicative identity; it equals zero in the zero ring."""
        return self.element(1)

    def hom(self, images, target):
        """Define a coefficient-preserving homomorphism by the coordinate images."""
        from .rings import RingHomomorphism

        return RingHomomorphism(self, target, images)

    def quotient(self, *generators, groebner_backend="native", groebner_limits=None):
        """Form a further quotient, retaining the ambient polynomial coordinates."""
        if len(generators) == 1 and isinstance(generators[0], (list, tuple)):
            generators = tuple(generators[0])
        expressions = tuple(self.element(g).as_expr() for g in generators)
        return CoordinateRing(
            self.ambient.ideal(
                *self.ideal.generators,
                *expressions,
                groebner_backend=groebner_backend,
                groebner_limits=groebner_limits,
            )
        )

    def base_change(self, domain):
        """Extend the coefficient field through its specified exact embedding.

        This is field extension, not specialization of transcendental parameters.
        """
        from .rings import _coefficient_embedding

        domain = exact_field(domain)
        _coefficient_embedding(self.domain, domain)
        ambient = PolynomialRing(self.variables, domain, self.ambient.order)
        return CoordinateRing(self.ideal.change_ring(ambient))

    @cached_property
    def krull_dimension(self):
        """The Krull dimension of S/I, with -1 for the zero ring."""
        return self.ideal.krull_dimension

    @cached_property
    def vector_space_dimension(self):
        """Finite algebra length, or infinity for positive-dimensional algebras."""
        if self.krull_dimension > 0:
            return sp.oo
        if not self.variables:
            return 0 if self.ideal.is_unit else 1
        return self.finite_algebra().dimension

    @cached_property
    def is_reduced(self):
        """The field-relative is reduced property of this coordinate ring."""
        from .properties import is_reduced

        return is_reduced(self)

    @cached_property
    def is_irreducible(self):
        """The field-relative is irreducible property of this coordinate ring."""
        from .properties import is_irreducible

        return is_irreducible(self)

    @cached_property
    def is_integral(self):
        """The field-relative is integral property of this coordinate ring."""
        from .properties import is_integral

        return is_integral(self)

    @cached_property
    def is_equidimensional(self):
        """The field-relative is equidimensional property of this coordinate ring."""
        from .properties import is_equidimensional

        return is_equidimensional(self)

    @cached_property
    def is_unmixed(self):
        """The field-relative is unmixed property of this coordinate ring."""
        from .properties import is_unmixed

        return is_unmixed(self)

    @cached_property
    def codimension(self):
        """The field-relative codimension property of this coordinate ring."""
        from .properties import codimension

        return codimension(self)

    def fraction_field(self):
        """Return the rational-function field of this integral coordinate ring."""
        from .rational_functions import FractionField

        return FractionField(self)

    def localize(self, *denominators):
        """Invert the supplied functions with a parent-aware localization."""
        from .rational_functions import Localization

        return Localization(self, denominators)

    def normalize(self, **limits):
        """Return the integral-closure coordinate ring of a reduced affine algebra."""
        from .normalization import normalize

        return normalize(self, **limits)

    @property
    def is_normal(self):
        """Whether this coordinate ring is normal."""
        from .normalization import is_normal

        return is_normal(self)

    def finite_algebra(self, *, max_dimension=None):
        """Convert a finite coordinate algebra to algroots, rejecting positive dimension.

        A scalar coefficient field uses a dummy coordinate constrained to zero
        because the finite-algebra backend requires a nonempty coordinate list.
        """
        if self.krull_dimension > 0:
            raise ValueError(
                "a positive-dimensional coordinate ring is not a finite quotient algebra"
            )
        variables = self.variables or (sp.Dummy("scalar_coordinate"),)
        generators = (
            self.ideal.generators
            if self.variables
            else (1 if self.ideal.is_unit else variables[0],)
        )
        return QuotientAlgebra.from_polynomials(
            generators,
            variables,
            domain=self.domain,
            max_dimension=max_dimension,
            groebner_backend=self.ideal.groebner_backend,
            groebner_limits=self.ideal.groebner_limits,
        )


__all__ = ["PolynomialRing", "Ideal", "CoordinateRing"]
