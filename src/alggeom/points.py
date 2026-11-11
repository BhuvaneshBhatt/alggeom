"""Complete finite-field enumeration and bounded-height rational point searches."""

from functools import reduce
from itertools import product
from math import gcd

import sympy as sp
from algroots.domains import field_element
from algroots.finite_fields import field_elements, field_order, is_finite_field

from ._algebra import integer


def _rational_pool(height, limit, coordinates):
    """Count reduced fractions before materializing the height-bounded search box."""
    if not coordinates:
        return (sp.S.Zero,)
    if 2 * height + 1 > limit:
        raise RuntimeError("point search exceeded max_candidates")
    phi = list(range(height + 1))
    for p in range(2, height + 1):
        if phi[p] == p:
            for multiple in range(p, height + 1, p):
                phi[multiple] -= phi[multiple] // p
    if (4 * sum(phi[1:]) - 1) ** coordinates > limit:
        raise RuntimeError("point search exceeded max_candidates")
    return tuple(
        sorted(
            {sp.Rational(a, b) for b in range(1, height + 1) for a in range(-height, height + 1)}
        )
    )


def point_search(scheme, height=None, *, max_candidates=1000000):
    """Return all finite-field points, or all QQ points in a specified height box.

    Affine rational height is max(|numerator|,denominator) in each reduced
    coordinate. Projective height is the maximum absolute entry of a primitive
    integer representative. A bounded empty QQ search is not a proof that no
    rational point exists. Exhausting the candidate budget raises without
    returning a prefix as a complete bounded enumeration.
    """
    from .schemes import AffineScheme, ProjectiveScheme

    if not isinstance(scheme, (AffineScheme, ProjectiveScheme)):
        raise TypeError("an affine or projective scheme is required")
    limit = integer(max_candidates, "max_candidates", 1)
    if scheme.is_empty:
        return ()
    finite = is_finite_field(scheme.domain)
    projective = isinstance(scheme, ProjectiveScheme)
    n = len(scheme.variables)
    if n == 0:
        return (scheme.point(()),)
    if projective and n == 1:
        return (scheme.point((1,)),)
    if finite:
        q = field_order(scheme.domain)
        count = (q**n - 1) // (q - 1) if projective else q**n
        if count > limit:
            raise RuntimeError("finite-field enumeration exceeded max_candidates")
        pool = tuple(scheme.domain.to_sympy(v) for v in field_elements(scheme.domain))
        if projective:
            candidates = (
                tuple([sp.S.Zero] * pivot + [sp.S.One] + list(tail))
                for pivot in range(n)
                for tail in product(pool, repeat=n - pivot - 1)
            )
        else:
            candidates = product(pool, repeat=n)
    else:
        if scheme.domain != sp.QQ:
            raise NotImplementedError("height-bounded point search currently supports QQ")
        if height is None:
            raise ValueError("a positive height bound is required for rational point search")
        height = integer(height, "height", 1)
        if projective:
            if n == 1:
                candidates = ((sp.S.One,),)
            else:
                # Primitive representatives have a unique positive first nonzero entry.
                if 2 * height + 1 > limit:
                    raise RuntimeError("projective height search exceeded max_candidates")
                count = (
                    sum(
                        int(sp.mobius(d)) * ((2 * (height // d) + 1) ** n - 1)
                        for d in range(1, height + 1)
                    )
                    // 2
                )
                if count > limit:
                    raise RuntimeError("projective height search exceeded max_candidates")
                candidates = (
                    tuple(values)
                    for pivot in range(n)
                    for first in range(1, height + 1)
                    for tail in product(range(-height, height + 1), repeat=n - pivot - 1)
                    if reduce(gcd, values := ([0] * pivot + [first] + list(tail))) == 1
                )
        else:
            candidates = product(_rational_pool(height, limit, n), repeat=n)
    result = []
    for values in candidates:
        substitution = dict(zip(scheme.variables, values, strict=True))
        if all(
            field_element(g.subs(substitution, simultaneous=True), scheme.domain)
            == scheme.domain.zero
            for g in scheme.ideal.generators
        ):
            result.append(scheme.point(values))
    return tuple(result)


__all__ = ["point_search"]
