"""Finite coefficient algebras and exact Galois orbits of polynomial factors."""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

import sympy as sp
from algroots.quotient import QuotientAlgebra
from algroots.subfields import field_normal_form
from sympy.polys.polyerrors import BasePolynomialError

from ._validation import (
    DecompositionError,
    DegreeLimitError,
    constant_field,
    parse_relation,
    validate_integer,
    validate_symbols,
)


def _monic(relation):
    x, z = relation.gens
    return sp.Poly(relation.as_expr(), x, domain=relation.domain.frac_field(z)).monic()


def _fresh(label, forbidden):
    names = {s.name for s in forbidden}
    while label in names:
        label = "_" + label
    return sp.Dummy(label)


def _equations(relation, degree, coefficients):
    """Universal remainder equations for a monic factor of the requested degree."""
    x, z = relation.gens
    base = relation.domain.frac_field(z)
    ring = base.poly_ring(*coefficients)
    generic = x**degree + sum((c * x ** (degree - i - 1) for i, c in enumerate(coefficients)))
    remainder = sp.Poly(_monic(relation).as_expr(), x, domain=ring).rem(
        sp.Poly(generic, x, domain=ring), auto=False
    )
    return tuple(remainder.all_coeffs())


def _evaluate_vector(polynomial, matrix, one):
    """Evaluate a scalar polynomial at a multiplication matrix by Horner iteration."""
    result = one.scalarmul(matrix.domain.zero)
    for coefficient in polynomial.all_coeffs():
        result = matrix * result + one.scalarmul(matrix.domain.convert(coefficient))
    return result


def _divide(relation, factor, ambient):
    """Divide the monic input by a factor using exact extension-field coefficients."""
    x = relation.gens[0]
    left = list(_monic(relation).all_coeffs())
    right = sp.Poly(factor, x).all_coeffs()
    quotient = []
    while len(left) >= len(right):
        lead = field_normal_form(ambient, left[0])
        quotient.append(lead)
        for i, c in enumerate(right):
            left[i] = field_normal_form(ambient, left[i] - lead * c)
        left.pop(0)
    expression = sum((c * x ** (len(quotient) - i - 1) for i, c in enumerate(quotient)))
    return (sp.expand(expression), tuple(field_normal_form(ambient, c) for c in left))


@dataclass(frozen=True)
class FactorOrbit:
    """One monic-factor orbit, in its explicitly constructed coefficient field."""

    relation: sp.Poly
    degree: int
    ambient: QuotientAlgebra
    coefficients: tuple[sp.Expr, ...]
    cofactor: sp.Expr

    @property
    def factor(self):
        """The monic factor represented by this orbit."""
        x = self.relation.gens[0]
        return x**self.degree + sum(
            (c * x ** (self.degree - i - 1) for i, c in enumerate(self.coefficients))
        )

    @property
    def extension_degree(self):
        """Degree of the factor coefficient field over constants(input)."""
        return self.ambient.dimension

    def verify(self):
        """Replay field embeddings, division identities and claimed orbit coverage."""
        try:
            x, z = self.relation.gens
            validate_integer(self.degree, "degree", 1)
            if self.degree > self.relation.degree(x) or len(self.coefficients) != self.degree:
                return False
            if (
                self.ambient.domain != self.relation.domain.frac_field(z)
                or not self.ambient.is_field
            ):
                return False
            validate_symbols((x, z, *self.ambient.variables))
            if (
                parse_relation(self.relation.as_expr(), (x, z), self.relation.domain)
                != self.relation
            ):
                return False
            if self.coefficients != tuple(
                field_normal_form(self.ambient, c) for c in self.coefficients
            ):
                return False
            quotient, remainder = _divide(self.relation, self.factor, self.ambient)
            return not any(remainder) and quotient == self.cofactor
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


@dataclass(frozen=True)
class FactorDegree:
    """An exhaustive, cyclic presentation of all factors of one fixed degree."""

    relation: sp.Poly
    degree: int
    quotient: QuotientAlgebra
    parameter: sp.Symbol
    weights: tuple[int, ...]
    defining_polynomial: sp.Poly
    coordinate_polynomials: tuple[sp.Poly, ...]
    component_polynomials: tuple[sp.Poly, ...]
    orbits: tuple[FactorOrbit, ...]

    @property
    def geometric_factor_count(self):
        """Number of factors of this degree over an algebraic closure."""
        return comb(self.relation.degree(self.relation.gens[0]), self.degree)

    def verify(self):
        """Replay universal equations, inverse maps and irreducible orbit coverage."""
        try:
            x, z = self.relation.gens
            A = self.quotient
            base = self.relation.domain.frac_field(z)
            validate_integer(self.degree, "degree", 1)
            if not 1 <= self.degree <= self.relation.degree(x) or len(A.variables) != self.degree:
                return False
            validate_symbols((x, z, self.parameter, *A.variables))
            expected = self.geometric_factor_count
            if A.domain != base or A.dimension != expected:
                return False
            rebuilt = QuotientAlgebra.from_polynomials(
                _equations(self.relation, self.degree, A.variables), A.variables, domain=base
            )
            if (
                A.groebner_basis != rebuilt.groebner_basis
                or A.standard_exponents != rebuilt.standard_exponents
            ):
                return False
            if len(self.weights) != self.degree or any(
                isinstance(w, bool) or not isinstance(w, int) for w in self.weights
            ):
                return False
            p = self.defining_polynomial
            if (
                p.gens != (self.parameter,)
                or p.domain != base
                or p.degree() != expected
                or (p.LC() != 1)
                or (p.gcd(p.diff()).degree() != 0)
            ):
                return False
            if len(self.coordinate_polynomials) != self.degree:
                return False
            if any(
                q.gens != p.gens or q.domain != base or q.degree() >= expected
                for q in self.coordinate_polynomials
            ):
                return False
            linear = sum((w * c for w, c in zip(self.weights, A.variables, strict=True)))
            matrix = A.multiplication_matrix(linear, representation="domain")
            one = A.coordinate_vector(1).to_DM(domain=base)
            if not _evaluate_vector(p, matrix, one).is_zero_matrix:
                return False
            if any(
                (
                    _evaluate_vector(q, matrix, one).to_Matrix() != A.coordinate_vector(c)
                    for q, c in zip(self.coordinate_polynomials, A.variables, strict=True)
                )
            ):
                return False
            images = tuple(q.as_expr() for q in self.coordinate_polynomials)
            if (
                not sp.Poly(
                    sum((w * q for w, q in zip(self.weights, images, strict=True)))
                    - self.parameter,
                    self.parameter,
                    domain=base,
                )
                .rem(p)
                .is_zero
            ):
                return False
            factors = self.component_polynomials
            if not factors or len(factors) != len(self.orbits):
                return False
            if any(
                f.gens != p.gens or f.domain != base or f.LC() != 1 or (not f.is_irreducible)
                for f in factors
            ):
                return False
            if sp.Poly(sp.prod(f.as_expr() for f in factors), self.parameter, domain=base) != p:
                return False
            for f, orbit in zip(factors, self.orbits, strict=True):
                C = orbit.ambient
                if (
                    orbit.relation != self.relation
                    or orbit.degree != self.degree
                    or C.variables != (self.parameter,)
                    or (C.domain != base)
                ):
                    return False
                expected_field = QuotientAlgebra.from_polynomials(
                    [f.as_expr()], (self.parameter,), domain=base
                )
                if (
                    C.groebner_basis != expected_field.groebner_basis
                    or C.standard_exponents != expected_field.standard_exponents
                ):
                    return False
                if (
                    orbit.coefficients != tuple(C.normal_form(q) for q in images)
                    or not orbit.verify()
                ):
                    return False
            return True
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


@dataclass(frozen=True)
class FactorIssue:
    """An uncompleted factor degree, never a negative decomposition proof."""

    degree: int
    reason: str


@dataclass(frozen=True)
class FactorDiscovery:
    """Factor extensions for requested degrees, grouped by Galois orbit."""

    relation: sp.Poly
    requested_degrees: tuple[int, ...]
    levels: tuple[FactorDegree, ...]
    issues: tuple[FactorIssue, ...]
    max_factor_dimension: int
    max_separating_attempts: int

    @property
    def orbits(self):
        """All discovered factor orbits in requested degree order."""
        return tuple(orbit for level in self.levels for orbit in level.orbits)

    @property
    def complete(self):
        """Coverage of requested factor degrees, including every conjugate orbit."""
        return (
            not self.issues
            and tuple(level.degree for level in self.levels) == self.requested_degrees
        )

    def verify(self):
        """Replay field embeddings, division identities and claimed orbit coverage."""
        try:
            x, z = self.relation.gens
            if (
                parse_relation(self.relation.as_expr(), (x, z), self.relation.domain)
                != self.relation
            ):
                return False
            n = self.relation.degree(x)
            if _degrees(self.requested_degrees, n) != self.requested_degrees:
                return False
            validate_integer(self.max_factor_dimension, "max_factor_dimension", 1)
            validate_integer(self.max_separating_attempts, "max_separating_attempts", 1)
            done = tuple(level.degree for level in self.levels)
            missing = tuple(i.degree for i in self.issues)
            if len(set((*done, *missing))) != len(done) + len(missing) or set(
                (*done, *missing)
            ) != set(self.requested_degrees):
                return False
            if done != tuple(d for d in self.requested_degrees if d in done) or missing != tuple(
                d for d in self.requested_degrees if d in missing
            ):
                return False
            if any(
                level.relation != self.relation
                or level.geometric_factor_count > self.max_factor_dimension
                or (
                    not any(
                        level.weights == tuple(j**i for i in range(level.degree))
                        for j in range(1, self.max_separating_attempts + 1)
                    )
                )
                or (not level.verify())
                for level in self.levels
            ):
                return False
            for issue in self.issues:
                validate_integer(issue.degree, "issue degree", 1)
                if issue.reason == "factor_dimension_limit":
                    if comb(n, issue.degree) <= self.max_factor_dimension:
                        return False
                elif issue.reason != "separator_limit":
                    return False
            return True
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


def _degrees(degrees, n):
    values = tuple(range(1, n + 1)) if degrees is None else tuple(degrees)
    for d in values:
        validate_integer(d, "factor degree", 1)
        if d > n:
            raise DecompositionError("factor degree exceeds the polynomial degree")
    if tuple(sorted(set(values))) != values:
        raise DecompositionError("factor degrees must be distinct and increasing")
    return values


def _level(relation, degree, max_dimension, max_attempts):
    """Split the universal factor algebra into its Galois orbit coefficient fields."""
    x, z = relation.gens
    base = relation.domain.frac_field(z)
    coefficients = tuple(_fresh(f"factor_coefficient_{i}", (x, z)) for i in range(degree))
    parameter = _fresh("factor_extension", (x, z, *coefficients))
    A = QuotientAlgebra.from_polynomials(
        _equations(relation, degree, coefficients),
        coefficients,
        domain=base,
        max_dimension=max_dimension,
    )
    # In characteristic zero the squarefree degree-n input has binomial(n,d) monic factors.
    dimension = comb(relation.degree(x), degree)
    if A.dimension != dimension:
        raise ArithmeticError("universal monic-factor algebra has an incorrect dimension")
    for attempt in range(1, max_attempts + 1):
        weights = tuple(attempt**i for i in range(degree))
        linear = sum((w * c for w, c in zip(weights, coefficients, strict=True)))
        p = A.minimal_polynomial(linear, parameter)
        if p.degree() == dimension and p.gcd(p.diff()).degree() == 0:
            break
    else:
        return None
    matrix = A.multiplication_matrix(linear, representation="domain")
    vector = A.coordinate_vector(1).to_DM(domain=base)
    powers = []
    for _ in range(dimension):
        powers.append(vector)
        vector = matrix * vector
    # The separating form gives a power basis; invert it to recover every factor coefficient.
    change = powers[0].hstack(*powers[1:])
    targets = tuple(A.coordinate_vector(c).to_DM(domain=base) for c in coefficients)
    coordinates = (change.inv() * targets[0].hstack(*targets[1:])).to_list()
    images = tuple(
        sp.Poly.from_list(
            [base.to_sympy(coordinates[j][i]) for j in range(dimension - 1, -1, -1)],
            parameter,
            domain=base,
        )
        for i in range(degree)
    )
    _, factorization = p.factor_list()
    if any((m != 1 for _, m in factorization)):
        raise ArithmeticError("cyclic factor-algebra polynomial did not split squarefreely")
    factors = tuple(
        sorted(
            (f.monic() for f, _ in factorization),
            key=lambda f: (f.degree(), sp.default_sort_key(f.as_expr())),
        )
    )
    if sp.Poly(sp.prod(f.as_expr() for f in factors), parameter, domain=base) != p:
        raise ArithmeticError("normalized extension factors did not reconstruct their polynomial")
    orbits = []
    for f in factors:
        B = QuotientAlgebra.from_polynomials([f.as_expr()], (parameter,), domain=base)
        values = tuple(B.normal_form(q.as_expr()) for q in images)
        factor = x**degree + sum((c * x ** (degree - i - 1) for i, c in enumerate(values)))
        cofactor, remainder = _divide(relation, factor, B)
        if any(remainder):
            raise ArithmeticError("discovered factor failed exact polynomial division")
        orbits.append(FactorOrbit(relation, degree, B, values, cofactor))
    return FactorDegree(relation, degree, A, parameter, weights, p, images, factors, tuple(orbits))


def discover_factors(
    relation,
    variables,
    *,
    degrees=None,
    constants=sp.QQ,
    max_degree=6,
    max_factor_dimension=32,
    max_separating_attempts=64,
):
    """Construct coefficient extensions and all monic factor orbits of requested degrees.

    Completeness is over constants(input), up to its Galois action, not numerical
    root selection. Limits bound input degree, factor-algebra dimension and the
    search for a separating form; they are not wall-clock limits.
    """
    x, z = tuple(variables)
    validate_symbols((x, z))
    constants = constant_field(constants)
    validate_integer(max_degree, "max_degree", 1)
    validate_integer(max_factor_dimension, "max_factor_dimension", 1)
    validate_integer(max_separating_attempts, "max_separating_attempts", 1)
    F = parse_relation(relation, (x, z), constants)
    n = F.degree(x)
    if n > max_degree:
        raise DegreeLimitError(f"output degree {n} exceeds max_degree={max_degree}")
    requested = _degrees(degrees, n)
    levels, issues = ([], [])
    for degree in requested:
        if comb(n, degree) > max_factor_dimension:
            issues.append(FactorIssue(degree, "factor_dimension_limit"))
            continue
        level = _level(F, degree, max_factor_dimension, max_separating_attempts)
        if level is None:
            issues.append(FactorIssue(degree, "separator_limit"))
        else:
            levels.append(level)
    return FactorDiscovery(
        F, requested, tuple(levels), tuple(issues), max_factor_dimension, max_separating_attempts
    )


__all__ = ["FactorOrbit", "FactorDegree", "FactorIssue", "FactorDiscovery", "discover_factors"]
