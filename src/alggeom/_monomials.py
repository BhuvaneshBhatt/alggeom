"""Exact combinatorics of monomial ideals."""

from functools import cache


def minimal_exponents(exponents):
    """Return the minimal monomial generators under componentwise divisibility."""
    values = tuple(sorted(set(tuple(e) for e in exponents)))
    return tuple(
        candidate
        for candidate in values
        if not any(
            divisor != candidate and all(a <= b for a, b in zip(divisor, candidate, strict=True))
            for divisor in values
        )
    )


def dimension(exponents, variable_count):
    """Compute dimension as n minus the minimum vertex cover of monomial supports.

    A coordinate prime contains a monomial precisely when it contains one of
    its variables. Branch on one support: every covering prime must choose at
    least one of those variables. Memoization shares the residual subproblems.
    """
    exponents = minimal_exponents(exponents)
    supports = frozenset(frozenset(i for i, power in enumerate(e) if power) for e in exponents)
    if frozenset() in supports:
        return -1

    @cache
    def height(edges):
        if not edges:
            return 0
        edge = min(edges, key=lambda e: (len(e), tuple(sorted(e))))
        return min(1 + height(frozenset(e for e in edges if v not in e)) for v in edge)

    return variable_count - height(supports)


def hilbert_numerator(exponents, variable_count):
    """Return coefficients of the monomial quotient's Hilbert numerator by lcm degree.

    Inclusion-exclusion over divisibility cones is combined on the lcm lattice,
    avoiding a separate term for every subset when many lcms coincide.
    """
    lattice = {(0,) * variable_count: 1}
    for generator in minimal_exponents(exponents):
        updates = {}
        for monomial, coefficient in lattice.items():
            lcm = tuple(max(a, b) for a, b in zip(monomial, generator, strict=True))
            updates[lcm] = updates.get(lcm, 0) - coefficient
        for monomial, coefficient in updates.items():
            lattice[monomial] = lattice.get(monomial, 0) + coefficient
        lattice = {m: c for m, c in lattice.items() if c}
    coefficients = {}
    for monomial, coefficient in lattice.items():
        degree = sum(monomial)
        coefficients[degree] = coefficients.get(degree, 0) + coefficient
    return {d: c for d, c in coefficients.items() if c}
