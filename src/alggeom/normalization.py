"""Affine integral closure by radical test ideals and finite endomorphism rings."""

from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from ._algebra import integer
from .ideals import CoordinateRing, PolynomialRing
from .properties import defining_ideal


def _test_ideal(ring, limits):
    """Derive the radical Jacobian test ideal independently of stored replay data."""
    from .maps import singular_locus

    return singular_locus(ring.ideal, assume_equidimensional=True).radical(**limits)


def _test_pair(ring, limits):
    """Use the radical Jacobian ideal as a Grauert-Remmert test ideal."""
    test = _test_ideal(ring, limits)
    denominator = next((g for g in test.generators if not ring.ideal.contains(g)), None)
    if denominator is None:
        raise NotImplementedError("the model is not generically smooth over its constants")
    multiples = ring.ideal + ring.ambient.ideal(*(denominator * g for g in test.generators))
    endomorphisms = multiples.quotient(test)
    return test, denominator, endomorphisms


def _domain_normal(ring, limits):
    _, g, endomorphisms = _test_pair(ring, limits)
    return endomorphisms == ring.ideal + ring.ambient.ideal(g)


def is_normal(value, **limits):
    """Decide normality over characteristic-zero fields, including reduced disconnected rings.

    Distinct irreducible components of a normal Noetherian ring are disjoint.
    Each component is then tested by equality A=End_A(J) for a radical
    Jacobian test ideal J. The zero ring is normal by the empty-scheme convention.
    """
    from .schemes import ProjectiveScheme

    if isinstance(value, ProjectiveScheme):
        return all(
            is_normal(value.affine_chart(i).scheme, **limits) for i in range(len(value.variables))
        )
    ideal = defining_ideal(value)
    if ideal.ring.domain.characteristic() != 0:
        raise NotImplementedError("normality currently requires characteristic-zero constants")
    if ideal.is_unit:
        return True
    if ideal.radical(**limits) != ideal:
        return False
    primes = ideal.minimal_primes(**limits)
    if any(not (p + q).is_unit for p, q in combinations(primes, 2)):
        return False
    return all(_domain_normal(CoordinateRing(p), limits) for p in primes)


@dataclass(frozen=True)
class NormalizationStep:
    """One finite birational extension End_A(J)=(gJ:J)/g, with replay data."""

    source: CoordinateRing
    target: CoordinateRing
    test_ideal: object
    denominator: sp.Expr
    numerators: tuple[sp.Expr, ...]

    def verify(self):
        """Replay the test ideal, colon module and saturated graph presentation."""
        if (
            self.source.domain != self.target.domain
            or self.target.variables[: len(self.source.variables)] != self.source.variables
        ):
            return False
        test = _test_ideal(self.source, {})
        if (
            test != self.test_ideal
            or self.source.ideal.contains(self.denominator)
            or not test.contains(self.denominator)
        ):
            return False
        g = self.denominator
        # A nonzero zero divisor can remove an entire component when the graph
        # is saturated. Birational division requires a regular denominator.
        if self.source.ideal.quotient(g) != self.source.ideal:
            return False
        expected = (
            self.source.ideal + self.source.ambient.ideal(*(g * j for j in test.generators))
        ).quotient(test)
        if expected != self.source.ideal + self.source.ambient.ideal(g, *self.numerators):
            return False
        variables = self.target.variables[len(self.source.variables) :]
        if len(variables) != len(self.numerators):
            return False
        graph = self.target.ambient.ideal(
            *self.source.ideal.generators,
            *(g * v - c for v, c in zip(variables, self.numerators, strict=True)),
        ).saturate(g)
        return graph == self.target.ideal


@dataclass(frozen=True)
class IntegralClosure:
    """A normalization ring and its finite birational coordinate-ring map.

    For a domain, inverse_functions express all normalized coordinates in the
    original function field, and module_generators span the integral closure
    over the input ring. For a reduced ring with several components, components
    contains their domain normalizations and ring is their explicit product.
    """

    mapping: object
    steps: tuple[NormalizationStep, ...] = ()
    inverse_functions: tuple = ()
    module_generators: tuple = ()
    components: tuple = ()

    @property
    def source(self):
        """The original reduced coordinate ring."""
        return self.mapping.source

    @property
    def ring(self):
        """The integral-closure coordinate ring."""
        return self.mapping.target

    @property
    def conductor(self):
        """The conductor ideal in the original ambient ring, for a domain input."""
        if self.components:
            raise NotImplementedError("the conductor of a product normalization is not provided")
        result = self.source.ambient.ideal(1)
        for element in self.module_generators:
            result = result.intersection(
                (self.source.ideal + self.source.ambient.ideal(element.denominator)).quotient(
                    element.numerator
                )
            )
        return result

    def verify(self):
        """Check finite-extension steps, birational coordinates and final normality."""
        if self.source.ideal.is_unit:
            return self.ring.ideal.is_unit
        if self.source.ideal.radical() != self.source.ideal:
            return False
        if self.components:
            if not all(c.verify() for c in self.components):
                return False
            primes = self.source.ideal.minimal_primes()
            if len(primes) != len(self.components) or not all(
                any(c.source.ideal == p for c in self.components) for p in primes
            ):
                return False
            expected = _product_closure(self.source, self.components)
            if len(expected.ring.variables) != len(self.ring.variables):
                return False
            rename = dict(zip(expected.ring.variables, self.ring.variables, strict=True))
            ideal = self.ring.ambient.ideal(
                *(g.xreplace(rename) for g in expected.ring.ideal.generators)
            )
            return ideal == self.ring.ideal and all(
                self.ring(f.as_expr().xreplace(rename)) == g
                for f, g in zip(expected.mapping.images, self.mapping.images, strict=True)
            )
        from .rational_functions import FractionField

        field = FractionField(self.source)
        if len(self.inverse_functions) != len(self.ring.variables):
            return False
        substitution = dict(
            zip(self.ring.variables, [f.as_expr() for f in self.inverse_functions], strict=True)
        )
        if any(
            field(g.subs(substitution, simultaneous=True)) != 0 for g in self.ring.ideal.generators
        ):
            return False
        if any(
            field(image.as_expr().subs(substitution, simultaneous=True)) != field(v)
            for image, v in zip(self.mapping.images, self.source.variables, strict=True)
        ):
            return False
        expected_module = (field.one,)
        current = self.source
        for step in self.steps:
            if step.source != current or not step.verify():
                return False
            fractions = self.inverse_functions[
                len(step.source.variables) : len(step.target.variables)
            ]
            expected_module = tuple(
                old * new for old in expected_module for new in (field.one, *fractions)
            )
            current = step.target
        return (
            current == self.ring
            and self.module_generators == expected_module
            and _domain_normal(self.ring, {})
        )


def _domain_closure(original, max_steps, max_new_variables, max_module_generators, limits):
    """Iterate endomorphism rings, transporting finite generators into the original field.

    Each saturated graph is birational, and products of the successive module
    generators span the full tower over the original coordinate ring.
    """
    from .rational_functions import FractionField

    field = FractionField(original)
    inverse = tuple(field(v) for v in original.variables)
    module = (field.one,)
    current = original
    steps = []
    for index in range(max_steps + 1):
        test, g, endomorphisms = _test_pair(current, limits)
        existing = current.ideal + current.ambient.ideal(g)
        if endomorphisms == existing:
            return IntegralClosure(
                original.hom(original.variables, current), tuple(steps), inverse, module
            )
        if index == max_steps:
            raise RuntimeError("integral closure exceeded max_steps")
        numerators = []
        span = existing
        for c in endomorphisms.groebner_basis:
            if not span.contains(c):
                numerators.append(c)
                span = span + current.ambient.ideal(c)
        new = tuple(sp.Dummy("integral_coordinate") for _ in numerators)
        if len(current.variables) + len(new) - len(original.variables) > max_new_variables:
            raise RuntimeError("integral closure exceeded max_new_variables")
        ambient = PolynomialRing((*current.variables, *new), current.domain)
        ideal = ambient.ideal(
            *current.ideal.generators, *(g * v - c for v, c in zip(new, numerators, strict=True))
        ).saturate(g)
        target = CoordinateRing(ideal)
        substitution = dict(zip(current.variables, [f.as_expr() for f in inverse], strict=True))
        fractions = tuple(field((c / g).subs(substitution, simultaneous=True)) for c in numerators)
        inverse = (*inverse, *fractions)
        # Products of finite module generators span through a tower of extensions.
        if len(module) * (len(fractions) + 1) > max_module_generators:
            raise RuntimeError("integral closure exceeded max_module_generators")
        module = tuple(old * new for old in module for new in (field.one, *fractions))
        steps.append(NormalizationStep(current, target, test, g, tuple(numerators)))
        current = target
    raise ArithmeticError("integral closure iteration did not terminate")


def _product_closure(original, components):
    """Present a product by orthogonal idempotents, retaining the original ring map."""
    idempotents = tuple(sp.Dummy("component") for _ in components)
    blocks = [tuple(sp.Dummy(v.name) for v in c.ring.variables) for c in components]
    ambient = PolynomialRing(
        (*idempotents, *(v for block in blocks for v in block)), original.domain
    )
    equations = [
        sum(idempotents) - 1,
        *(e * e - e for e in idempotents),
        *(e * f for e, f in combinations(idempotents, 2)),
    ]
    images = [sp.S.Zero] * len(original.variables)
    for i, component in enumerate(components):
        e = idempotents[i]
        rename = dict(zip(component.ring.variables, blocks[i], strict=True))
        equations.extend((1 - e) * v for v in blocks[i])
        equations.extend(e * g.xreplace(rename) for g in component.ring.ideal.generators)
        for j, image in enumerate(component.mapping.images):
            images[j] += e * image.as_expr().xreplace(rename)
    ring = CoordinateRing(ambient.ideal(*equations))
    return IntegralClosure(original.hom(images, ring), components=tuple(components))


def integral_closure(
    value,
    *,
    max_steps=16,
    max_new_variables=64,
    max_module_generators=1024,
    max_components=32,
    **decomposition_limits,
):
    """Compute the affine normalization with a finite birational map and exact evidence.

    Inputs must be reduced algebras over supported characteristic-zero fields.
    Domains use the Grauert-Remmert endomorphism-ring algorithm. Several
    components are normalized separately and assembled as a product. Limits
    raise without returning an intermediate ring as a completed normalization.
    """
    from .schemes import ProjectiveScheme

    if isinstance(value, ProjectiveScheme):
        raise NotImplementedError(
            "normalize affine charts explicitly; a projective normalized model is not constructed"
        )
    ideal = defining_ideal(value)
    if ideal.ring.domain.characteristic() != 0:
        raise NotImplementedError(
            "integral closure currently requires characteristic-zero constants"
        )
    steps = integer(max_steps, "max_steps")
    variables = integer(max_new_variables, "max_new_variables", 1)
    generators = integer(max_module_generators, "max_module_generators", 1)
    if ideal.radical(**decomposition_limits) != ideal:
        raise ValueError("normalization requires a reduced ring; reduce explicitly first")
    original = CoordinateRing(ideal)
    if ideal.is_unit:
        return IntegralClosure(original.hom(original.gens, original))
    primes = ideal.minimal_primes(**decomposition_limits)
    if len(primes) > integer(max_components, "max_components", 1):
        raise RuntimeError("integral closure exceeded max_components")
    if len(primes) == 1:
        return _domain_closure(original, steps, variables, generators, decomposition_limits)
    components = [
        _domain_closure(CoordinateRing(p), steps, variables, generators, decomposition_limits)
        for p in primes
    ]
    return _product_closure(original, components)


def normalize(value, **limits):
    """Return the integral-closure coordinate ring of a reduced affine algebra."""
    return integral_closure(value, **limits).ring


def normalization_map(value, **limits):
    """Return the normalization morphism for an affine scheme, or its ring homomorphism."""
    from .schemes import AffineScheme

    data = integral_closure(value, **limits)
    if isinstance(value, AffineScheme):
        return AffineScheme(data.ring).morphism(
            value, tuple(x.as_expr() for x in data.mapping.images)
        )
    return data.mapping


__all__ = [
    "NormalizationStep",
    "IntegralClosure",
    "is_normal",
    "integral_closure",
    "normalize",
    "normalization_map",
]
