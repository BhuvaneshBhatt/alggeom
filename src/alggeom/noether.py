"""Constructive Noether normalization with finite-module and integrality data."""

from dataclasses import dataclass
from itertools import islice, product

import sympy as sp
from algroots.quotient import standard_exponent_count

from ._algebra import integer
from .ideals import CoordinateRing, PolynomialRing
from .properties import defining_ideal


def _finite_data(ring, dimension, max_generators):
    """Read triangular monic relations and finite spanning monomials from lex order."""
    count = len(ring.variables) - dimension
    bounds = []
    relations = []
    for i in range(count):
        candidates = [g for g in ring.ideal._basis if g.LM[i] > 0 and sum(g.LM) == g.LM[i]]
        if not candidates:
            raise ArithmeticError("the coordinate change has not certified a finite projection")
        relation = min(candidates, key=lambda g: g.LM[i])
        bounds.append(relation.LM[i])
        relations.append(relation.as_expr())
    leading = [g.LM[:count] for g in ring.ideal._basis if not any(g.LM[count:])]
    monomials = []
    for powers in product(*(range(b) for b in bounds)):
        if any(all(a <= b for a, b in zip(lm, powers, strict=True)) for lm in leading):
            continue
        monomials.append(sp.prod(v**e for v, e in zip(ring.variables[:count], powers, strict=True)))
        if len(monomials) > max_generators:
            raise RuntimeError("Noether module generators exceeded max_generators")
    return tuple(relations), tuple(monomials)


def _monic_direction(polynomial, variables, domain, attempts, max_degree):
    """Try linear projections, then a degree-separated triangular automorphism.

    Distinct base-B exponent weights make every monomial of the selected
    equation have a different highest power of the eliminated coordinate.
    This fallback works over finite fields where linear projections may fail.
    """
    variable = variables[0]
    degree = sp.Poly(polynomial, *variables).total_degree()
    characteristic = int(domain.characteristic())
    if characteristic:
        from algroots.finite_fields import field_elements

        candidates = tuple(domain.to_sympy(c) for c in islice(field_elements(domain), attempts))
    else:
        candidates = range(-2, 3)
    for attempt, coefficients in enumerate(product(candidates, repeat=len(variables) - 1)):
        if attempt >= attempts:
            break
        substitution = {
            v: v + c * variable for v, c in zip(variables[1:], coefficients, strict=True)
        }
        converted = PolynomialRing(variables, domain).convert(
            polynomial.subs(substitution, simultaneous=True)
        )
        transformed = sp.Poly(converted, variable)
        if (
            not (transformed.LC().free_symbols & set(variables))
            and transformed.degree() <= max_degree
        ):
            return substitution
    weight = degree + 1
    powers = tuple(weight**i for i in range(1, len(variables)))
    if degree * max((1, *powers)) > max_degree:
        raise RuntimeError("Noether coordinate degree exceeded max_transform_degree")
    return {v: v + variable**power for v, power in zip(variables[1:], powers, strict=True)}


@dataclass(frozen=True)
class NoetherNormalization:
    """A finite polynomial subring and an explicit isomorphic coordinate presentation.

    coordinate_map goes from the input ring to the transformed ring. The
    polynomial base embeds through base_embedding. Integral relations are
    triangular and monic in the nonbase coordinates. module_generators span
    over the polynomial base; they are not asserted to be a free basis.
    """

    coordinate_map: object
    inverse_map: object
    base_embedding: object
    integral_relations: tuple[sp.Expr, ...]
    module_generators: tuple[sp.Expr, ...]

    @property
    def ring(self):
        """The transformed coordinate ring, isomorphic to the input."""
        return self.coordinate_map.target

    @property
    def base(self):
        """The embedded polynomial coordinate ring."""
        return self.base_embedding.source

    @property
    def projection(self):
        """The finite affine-scheme morphism to the polynomial base."""
        from .schemes import AffineScheme

        return AffineScheme(self.ring).morphism(AffineScheme(self.base), self.base.variables)

    @property
    def generic_rank(self):
        """The finite vector-space length after extending the polynomial base to its field."""
        variables = self.ring.variables[: len(self.ring.variables) - len(self.base.variables)]
        if not variables:
            return 1
        # Nonbase coordinates precede base coordinates in this lex basis.
        # Localizing the base makes its nonzero coefficients invertible, so
        # projecting each leading monomial gives the generic leading ideal.
        # This reuses the verified basis even in positive characteristic.
        count = len(variables)
        leading = tuple(g.LM[:count] for g in self.ring.ideal._basis)
        return standard_exponent_count(leading, count)

    def verify(self):
        """Check inverse coordinates, base injection and the monic finite-module certificate."""
        if (
            self.inverse_map.source != self.ring
            or self.inverse_map.target != self.coordinate_map.source
            or self.base_embedding.target != self.ring
            or not self.base.ideal.is_zero
        ):
            return False
        count = len(self.ring.variables) - len(self.base.variables)
        if (
            count < 0
            or self.base.variables != self.ring.variables[count:]
            or self.base_embedding.images != tuple(self.ring(v) for v in self.base.variables)
        ):
            return False
        if self.coordinate_map.compose(self.inverse_map).images != self.ring.gens:
            return False
        if self.inverse_map.compose(self.coordinate_map).images != self.coordinate_map.source.gens:
            return False
        if not self.base_embedding.kernel.is_zero:
            return False
        if len(self.integral_relations) != count:
            return False
        for i, equation in enumerate(self.integral_relations):
            if not self.ring.ideal.contains(equation) or equation.free_symbols & set(
                self.ring.variables[:i]
            ):
                return False
            polynomial = sp.Poly(equation, self.ring.variables[i])
            if polynomial.degree() < 1 or polynomial.LC() != 1:
                return False
        if not self.module_generators:
            return False
        try:
            _, monomials = _finite_data(
                self.ring, len(self.base.variables), len(self.module_generators)
            )
        except RuntimeError:
            return False
        return monomials == self.module_generators


def noether_normalization(
    value, *, names=None, max_linear_attempts=64, max_transform_degree=256, max_generators=1024
):
    """Construct a certified finite projection to affine space over the coefficient field.

    Successive coordinate automorphisms make an elimination equation monic.
    The fallback uses triangular polynomial changes when the field is too small
    for a suitable linear direction. Resource exhaustion raises without
    claiming that an uncertified generically finite projection is finite.
    """
    ideal = defining_ideal(value)
    if ideal.is_unit:
        raise ValueError("Noether normalization requires a nonzero coordinate ring")
    attempts = integer(max_linear_attempts, "max_linear_attempts", 1)
    limit = integer(max_transform_degree, "max_transform_degree", 1)
    maximum = integer(max_generators, "max_generators", 1)
    original = CoordinateRing(ideal)
    variables = (
        tuple(sp.Dummy("normalized_" + v.name) for v in ideal.ring.variables)
        if names is None
        else tuple(names)
    )
    if len(variables) != len(ideal.ring.variables):
        raise ValueError("one normalization coordinate is required per ambient variable")
    ambient = PolynomialRing(variables, ideal.ring.domain, "lex")
    rename = dict(zip(ideal.ring.variables, variables, strict=True))
    current = ambient.ideal(*(g.xreplace(rename) for g in ideal.generators))
    images = list(variables)
    inverse = list(ideal.ring.variables)
    dimension = ideal.krull_dimension
    for index in range(len(variables) - dimension):
        eliminated = current.eliminate(*variables[:index]) if index else current
        equation = min(
            eliminated.generators, key=lambda f: sp.Poly(f, *variables[index:]).total_degree()
        )
        substitution = _monic_direction(
            equation, variables[index:], ambient.domain, attempts, limit
        )
        old_inverse = dict(zip(variables, inverse, strict=True))
        inverse_change = {v: 2 * v - substitution[v] for v in substitution}
        inverse = [inverse_change.get(v, v).xreplace(old_inverse) for v in variables]
        images = [sp.expand(f.subs(substitution, simultaneous=True)) for f in images]
        current = ambient.ideal(
            *(sp.expand(g.subs(substitution, simultaneous=True)) for g in current.generators)
        )
    normalized = CoordinateRing(current)
    forward = original.hom(images, normalized)
    backward = normalized.hom(inverse, original)
    base = PolynomialRing(
        variables[len(variables) - dimension :] if dimension else (), ambient.domain
    ).quotient()
    embedding = base.hom(base.variables, normalized)
    relations, monomials = _finite_data(normalized, dimension, maximum)
    result = NoetherNormalization(forward, backward, embedding, relations, monomials)
    if not result.verify():
        raise ArithmeticError("Noether normalization certificate did not verify")
    return result


__all__ = ["NoetherNormalization", "noether_normalization"]
