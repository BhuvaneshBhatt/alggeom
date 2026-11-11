"""Exact radicals and primary components with coefficient-field-specific algorithms."""

from dataclasses import dataclass
from itertools import combinations

import sympy as sp
from algroots.domains import polynomial
from algroots.finite_decomposition import primary_decomposition as finite_decomposition
from algroots.finite_decomposition import reduced_algebra
from algroots.quotient import QuotientAlgebra, standard_exponent_count

from ._algebra import integer
from .ideals import CoordinateRing, Ideal, PolynomialRing


class DecompositionLimitError(RuntimeError):
    """An ideal decomposition exceeds a declared resource limit; no partial answer is returned."""


@dataclass
class _Budget:
    dimension: int
    nodes: int
    exponent: int
    attempts: int
    used: int = 0

    def visit(self):
        self.used += 1
        if self.used > self.nodes:
            raise DecompositionLimitError("ideal decomposition exceeded max_nodes")


def _budget(max_dimension, max_nodes, max_exponent, max_separating_attempts):
    return _Budget(
        *(
            integer(v, n, 1)
            for v, n in zip(
                (max_dimension, max_nodes, max_exponent, max_separating_attempts),
                ("max_dimension", "max_nodes", "max_exponent", "max_separating_attempts"),
                strict=True,
            )
        )
    )


def _intersection(ring, ideals):
    result = ring.ideal(1)
    for ideal in ideals:
        result = result.intersection(ideal)
    return result


def _finite_pairs(ideal, budget, primary):
    """Reduce dimension-zero arithmetic to the finite-algebra backend."""
    size = standard_exponent_count([g.LM for g in ideal._basis], len(ideal.ring.variables))
    if size > budget.dimension:
        raise DecompositionLimitError("generic finite-algebra dimension exceeded max_dimension")
    algebra = QuotientAlgebra.from_polynomials(
        ideal.groebner_basis,
        ideal.ring.variables,
        domain=ideal.ring.domain,
        max_dimension=budget.dimension,
        groebner_backend=ideal.groebner_backend,
        groebner_limits=ideal.groebner_limits,
    )
    options = dict(groebner_backend=ideal.groebner_backend, groebner_limits=ideal.groebner_limits)
    if not primary:
        algebra = reduced_algebra(algebra, **options)
    pairs = finite_decomposition(algebra, max_separating_attempts=budget.attempts, **options)
    return [
        (
            ideal._generated_ideal(p.as_expr() for p in component.groebner_basis.polys),
            ideal._generated_ideal(p.as_expr() for p in residue.groebner_basis.polys),
        )
        for component, residue in pairs
    ]


def _independent_coordinates(ideal):
    supports = tuple(frozenset(i for i, p in enumerate(g.LM) if p) for g in ideal._basis)
    dimension = ideal.krull_dimension
    for indices in combinations(range(len(ideal.ring.variables)), dimension):
        if not any(support <= frozenset(indices) for support in supports):
            return tuple(ideal.ring.variables[i] for i in indices)
    raise RuntimeError("the initial monomial ideal has no dimension-sized independent set")


def _generic_algebra(ideal, independent):
    """Localize a maximal independent coordinate set and record contraction denominators.

    Lex order places dependent coordinates before independent coordinates.
    The product h of their leading coefficients is sufficient to contract
    the generic ideal: I*K(U)[V] intersect K[U,V] = I:h^infinity.
    Existing coefficient-field parameters remain coefficient parameters.
    """
    dependent = tuple(v for v in ideal.ring.variables if v not in independent)
    field = ideal.ring.domain
    ground = field.domain if field.is_FractionField else field
    parameters = field.symbols if field.is_FractionField else ()
    generic_field = ground.frac_field(*parameters, *independent)
    block = PolynomialRing((*dependent, *independent), field, "lex")
    basis = ideal.change_ring(block).groebner_basis
    leading = tuple(polynomial(g, dependent, generic_field).LC() for g in basis)
    h = ideal.ring.convert(sp.prod(leading))
    generic = PolynomialRing(dependent, generic_field, "lex").ideal(*basis)
    return generic, h


def _contract(ideal, ambient):
    """Contract a monic generic Groebner basis by clearing and saturating denominators."""
    numerators = []
    denominators = []
    for generator in ideal.groebner_basis:
        numerator, denominator = sp.together(generator).as_numer_denom()
        numerators.append(ambient.convert(numerator))
        denominators.append(ambient.convert(denominator))
    return ambient.ideal(*numerators).saturate(sp.prod(denominators))


def _saturation_exponent(ideal, h, saturated, limit):
    """Find N with h^N*(I:h^infinity) contained in I, certifying the primary split."""
    if saturated == ideal:
        return 0
    exponent = 1
    while exponent <= limit:
        if all(ideal.contains(h**exponent * g) for g in saturated.generators):
            return exponent
        exponent = (
            2 * exponent if 2 * exponent <= limit else limit if exponent < limit else limit + 1
        )
    raise DecompositionLimitError("saturation stabilization exceeded max_exponent")


def _merge(pairs, ring, primary):
    groups = []
    for component, prime in pairs:
        for i, (old, associated) in enumerate(groups):
            if prime == associated:
                groups[i] = (old.intersection(component), prime)
                break
        else:
            groups.append((component, prime))
    if not primary:
        return [(p, p) for _, p in groups if not any(q != p and q.is_subset(p) for _, q in groups)]
    index = 0
    while index < len(groups):
        others = _intersection(ring, [q for j, (q, _) in enumerate(groups) if j != index])
        if others.is_subset(groups[index][0]):
            groups.pop(index)
        else:
            index += 1
    return groups


def _decompose(ideal, budget, primary, memo):
    key = tuple(ideal.groebner_basis)
    if key in memo:
        return memo[key]
    budget.visit()
    if ideal.is_unit:
        result = []
    elif ideal.is_zero:
        result = [(ideal, ideal)]
    elif len(ideal._basis) == 1:
        # In a polynomial UFD, powers of irreducible principal ideals are primary.
        generator = polynomial(ideal.groebner_basis[0], ideal.ring.variables, ideal.ring.domain)
        result = [
            (ideal.ring.ideal(f.as_expr() ** (e if primary else 1)), ideal.ring.ideal(f.as_expr()))
            for f, e in generator.factor_list()[1]
        ]
    elif ideal.krull_dimension == 0:
        result = _finite_pairs(ideal, budget, primary)
    else:
        independent = _independent_coordinates(ideal)
        generic, h = _generic_algebra(ideal, independent)
        finite = _finite_pairs(generic, budget, primary)
        result = [(_contract(q, ideal.ring), _contract(p, ideal.ring)) for q, p in finite]
        if primary:
            saturated = ideal.saturate(h)
            exponent = _saturation_exponent(ideal, h, saturated, budget.exponent)
        else:
            exponent = 1
        vertical = ideal + ideal.ring.ideal(h**exponent)
        if not vertical.is_unit:
            result.extend(_decompose(vertical, budget, primary, memo))
        result = _merge(result, ideal.ring, primary)
    memo[key] = result
    return result


def _pairs(ideal, primary, max_dimension, max_nodes, max_exponent, max_separating_attempts):
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    budget = _budget(max_dimension, max_nodes, max_exponent, max_separating_attempts)
    if ideal.ring.domain.characteristic() and ideal.ring.variables:
        if ideal.is_unit:
            return []
        if ideal.is_zero:
            return [(ideal, ideal)]
        if ideal.krull_dimension != 0:
            raise NotImplementedError(
                "positive-dimensional primary components currently require characteristic zero"
            )
        return _finite_pairs(ideal, budget, primary)

    if not ideal.ring.variables:
        return [] if ideal.is_unit else [(ideal, ideal)]
    return _decompose(ideal, budget, primary, {})


def primary_decomposition(
    ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64
):
    """Return an irredundant primary decomposition, including embedded components.

    Generic localization at independent coordinates reduces to finite algebras.
    Saturation contracts the generic primary factors; a certified exponent
    splits off the remaining closed locus. Components are over the declared
    field, not automatically over its algebraic closure. Resource exhaustion
    raises and never returns an apparently complete prefix. Decompositions of
    embedded components need not be unique.
    """
    return tuple(
        q
        for q, _ in _pairs(
            ideal, True, max_dimension, max_nodes, max_exponent, max_separating_attempts
        )
    )


def associated_primes(
    ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64
):
    """Return all associated primes, including embedded primes, over the declared field."""
    return tuple(
        p
        for _, p in _pairs(
            ideal, True, max_dimension, max_nodes, max_exponent, max_separating_attempts
        )
    )


def minimal_primes(
    ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64
):
    """Return the prime ideals of the reduced irreducible components over the declared field."""
    pairs = _pairs(ideal, False, max_dimension, max_nodes, max_exponent, max_separating_attempts)
    return tuple(p for _, p in _merge(pairs, ideal.ring, False))


def radical(
    ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64
):
    """Compute the radical by characteristic-zero components or finite-field Frobenius roots."""
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    if ideal.ring.domain.characteristic():
        from algroots.finite_fields import field_order, is_finite_field

        if not is_finite_field(ideal.ring.domain):
            raise NotImplementedError("positive-characteristic radicals require finite constants")

        from .rings import RingHomomorphism

        q = field_order(ideal.ring.domain)
        budget = _budget(max_dimension, max_nodes, max_exponent, max_separating_attempts)
        current = ideal
        for _ in range(budget.nodes):
            if current.is_zero or current.is_unit:
                return current
            quotient = CoordinateRing(current)
            source = ideal.ring.quotient()
            root = RingHomomorphism(
                source, quotient, tuple(v**q for v in ideal.ring.variables)
            ).kernel
            if root == current:
                return current
            current = root
        raise DecompositionLimitError("Frobenius radical exceeded max_nodes")
    primes = minimal_primes(
        ideal,
        max_dimension=max_dimension,
        max_nodes=max_nodes,
        max_exponent=max_exponent,
        max_separating_attempts=max_separating_attempts,
    )
    return _intersection(ideal.ring, primes)


def is_prime(ideal, **limits):
    """Decide whether a proper ideal is prime over its declared coefficient field."""
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    if ideal.is_unit:
        return False
    primes = minimal_primes(ideal, **limits)
    return len(primes) == 1 and ideal == primes[0]


def is_primary(ideal, **limits):
    """Decide whether a proper ideal is primary, retaining nilpotents."""
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    return not ideal.is_unit and len(primary_decomposition(ideal, **limits)) == 1


__all__ = [
    "DecompositionLimitError",
    "radical",
    "minimal_primes",
    "associated_primes",
    "primary_decomposition",
    "is_prime",
    "is_primary",
]
