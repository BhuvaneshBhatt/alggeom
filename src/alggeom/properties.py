"""Geometric predicates with explicit field-relative and empty-scheme conventions."""

import sympy as sp

from .ideals import CoordinateRing, Ideal, PolynomialRing


def defining_ideal(value):
    """Extract the defining ideal of a ring or embedded scheme."""
    if isinstance(value, Ideal):
        return value
    if isinstance(value, PolynomialRing):
        return value.ideal()
    if isinstance(value, CoordinateRing):
        return value.ideal
    from .schemes import AffineScheme, ProjectiveScheme

    if isinstance(value, (AffineScheme, ProjectiveScheme)):
        return value.ideal
    raise TypeError("an ideal, polynomial ring, coordinate ring or embedded scheme is required")


def is_reduced(value, **limits):
    """Decide whether the scheme or coordinate ring has no nonzero nilpotents."""
    ideal = defining_ideal(value)
    return ideal == ideal.radical(**limits)


def is_irreducible(value, **limits):
    """Decide irreducibility over the declared field; the empty scheme is not irreducible."""
    return len(defining_ideal(value).minimal_primes(**limits)) == 1


def is_integral(value, **limits):
    """Decide whether the coordinate ring is a nonzero domain over its declared field."""
    from .components import is_prime

    return is_prime(defining_ideal(value), **limits)


def is_equidimensional(value, **limits):
    """Decide equality of minimal-component dimensions; the empty scheme satisfies it."""
    return len({p.krull_dimension for p in defining_ideal(value).minimal_primes(**limits)}) <= 1


def is_unmixed(value, **limits):
    """Decide equality of all associated-prime dimensions, including embedded primes."""
    return len({p.krull_dimension for p in defining_ideal(value).associated_primes(**limits)}) <= 1


def codimension(value):
    """Return embedded codimension, with infinity for the empty scheme."""
    ideal = defining_ideal(value)
    return sp.oo if ideal.is_unit else len(ideal.ring.variables) - ideal.krull_dimension


def is_smooth(value, point=None):
    """Decide smoothness globally, or at a supplied exact field-valued point.

    Pointwise computations use the local tangent cone and embedding dimension.
    Global mixed-dimensional calculations use the reduced singular locus.
    """
    from .local import embedding_dimension, local_dimension
    from .maps import singular_locus

    if point is not None:
        return embedding_dimension(value, point) == local_dimension(value, point)
    result = singular_locus(value, reduced=True)
    return defining_ideal(result).is_unit


__all__ = [
    "is_reduced",
    "is_irreducible",
    "is_integral",
    "is_equidimensional",
    "is_unmixed",
    "codimension",
    "is_smooth",
]
