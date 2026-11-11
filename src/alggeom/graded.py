"""Standard-graded Hilbert invariants, dimensions and projective closure."""

import sympy as sp

from ._algebra import integer
from ._monomials import hilbert_numerator
from .ideals import CoordinateRing, Ideal, PolynomialRing


def _ideal(value):
    if isinstance(value, CoordinateRing):
        return value.ideal
    if isinstance(value, PolynomialRing):
        return value.ideal()
    if not isinstance(value, Ideal):
        raise TypeError("an ideal, coordinate ring, or polynomial ring is required")
    return value


def _parameter(ideal, parameter, name):
    forbidden = {
        *ideal.ring.variables,
        *(ideal.ring.domain.symbols if ideal.ring.domain.is_FractionField else ()),
    }
    if parameter is None:
        candidate = sp.Symbol(name)
        parameter = sp.Dummy(name) if candidate in forbidden else candidate
    if not isinstance(parameter, sp.Symbol) or parameter in forbidden:
        raise ValueError("a fresh Symbol is required for the Hilbert parameter")
    return parameter


def is_homogeneous(ideal):
    """Decide homogeneity for the standard grading, treating coefficient-field elements as scalars."""
    ideal = _ideal(ideal)
    return all(len({sum(m) for m in p.monoms()}) <= 1 for p in ideal._basis)


def _numerator(ideal):
    if not is_homogeneous(ideal):
        raise ValueError("Hilbert invariants require a homogeneous ideal in the standard grading")
    return hilbert_numerator((p.LM for p in ideal._basis), len(ideal.ring.variables))


def hilbert_series(ideal, *, parameter=None):
    """Return the exact rational Hilbert series of S/I, not of I as an S-module.

    The standard grading assigns degree one to each coordinate. The numerator
    is computed on the leading monomial ideal by lcm inclusion-exclusion.
    Arbitrary weights and multigradings are not inferred.
    """
    ideal = _ideal(ideal)
    parameter = _parameter(ideal, parameter, "t")
    coefficients = _numerator(ideal)
    numerator = sum(c * parameter**d for d, c in coefficients.items())
    return sp.cancel(numerator / (1 - parameter) ** len(ideal.ring.variables))


def hilbert_function(ideal, degree):
    """Return dim_K(S/I)_degree for a nonnegative integer degree."""
    ideal = _ideal(ideal)
    degree = integer(degree, "degree")
    coefficients = _numerator(ideal)
    count = len(ideal.ring.variables)
    if not count:
        return coefficients.get(degree, 0)
    return int(
        sum(
            c * sp.binomial(degree - d + count - 1, count - 1)
            for d, c in coefficients.items()
            if d <= degree
        )
    )


def hilbert_polynomial(ideal, *, parameter=None):
    """Return the eventual polynomial of the graded-piece dimensions as a QQ Poly.

    This is not the cumulative Hilbert-Samuel polynomial. An Artinian graded
    quotient has a finite Hilbert series and eventual Hilbert polynomial zero.
    """
    ideal = _ideal(ideal)
    parameter = _parameter(ideal, parameter, "n")
    coefficients = _numerator(ideal)
    count = len(ideal.ring.variables)
    expression = (
        0
        if not count
        else sum(
            c * sp.prod(parameter - d + j for j in range(1, count)) / sp.factorial(count - 1)
            for d, c in coefficients.items()
        )
    )
    return sp.Poly(sp.expand(expression), parameter, domain=sp.QQ)


def krull_dimension(ideal):
    """Return dim(S/I), also for inhomogeneous ideals; -1 denotes the zero ring."""
    return _ideal(ideal).krull_dimension


def projective_dimension(ideal):
    """Return dim Proj(S/I), with -1 for an empty projective scheme."""
    ideal = _ideal(ideal)
    _numerator(ideal)
    return max(-1, ideal.krull_dimension - 1)


def degree(ideal):
    """Return standard-graded Hilbert multiplicity; an Artinian quotient gives its length.

    For a nonempty projective scheme this is its degree. The zero ring gives
    zero. No smoothness, radicality or equidimensionality is assumed.
    """
    ideal = _ideal(ideal)
    coefficients = _numerator(ideal)
    dimension = ideal.krull_dimension
    if dimension < 0:
        return 0
    q = sp.Dummy("hilbert_degree")
    numerator = sum(c * q**d for d, c in coefficients.items())
    codimension = len(ideal.ring.variables) - dimension
    return int(
        (-1) ** codimension
        * sp.diff(numerator, q, codimension).subs(q, 1)
        / sp.factorial(codimension)
    )


def projective_closure(ideal, *, homogenizer=None):
    """Return the homogenized ideal saturated at a new projective coordinate.

    Saturation removes components at infinity introduced by homogenizing only
    the supplied affine generators. The resulting ring appends the homogenizer
    to the original coordinates. Scheme structure is retained.
    """
    ideal = _ideal(ideal)
    if homogenizer is None:
        homogenizer = sp.Dummy("homogenizer")
    if not isinstance(homogenizer, sp.Symbol) or homogenizer in ideal.ring.variables:
        raise ValueError("the homogenizer must be a fresh Symbol")
    ring = PolynomialRing((*ideal.ring.variables, homogenizer), ideal.ring.domain, ideal.ring.order)
    generators = []
    for g in ideal.generators:
        poly = ideal.ring._poly(g)
        total = max(map(sum, poly.monoms()))
        generators.append(
            sum(
                ideal.ring.domain.to_sympy(c)
                * sp.prod(v**p for v, p in zip(ideal.ring.variables, m, strict=True))
                * homogenizer ** (total - sum(m))
                for m, c in poly.items()
            )
        )
    return ring.ideal(*generators).saturate(homogenizer)


__all__ = [
    "is_homogeneous",
    "hilbert_series",
    "hilbert_function",
    "hilbert_polynomial",
    "krull_dimension",
    "projective_dimension",
    "degree",
    "projective_closure",
]
