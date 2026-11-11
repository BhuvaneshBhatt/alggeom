"""Recursive Gröbner covers over prime parameter ideals.

A branch works in Frac(QQ[a]/P) for a certified prime ideal P. Exceptional
leading coefficients enlarge P; their minimal primes are recursively visited.
The covers describe algebraic-closure parameter points, but the coefficient
fraction field and associated primes are computed over QQ.
"""

from dataclasses import dataclass

import sympy as sp

from .cover_proofs import check_parameter_partition
from .ideals import PolynomialRing
from .parameter_fields import (
    QuotientGroebnerResult,
    groebner_over_prime_parameters,
)


class ParameterCoverLimitError(RuntimeError):
    """A recursive algebraic parameter cover exceeds the declared budget."""


@dataclass(frozen=True)
class PrimeCoverPiece:
    """Constructible locus and a generic-basis proof valid everywhere on it."""

    condition: sp.Expr
    certificate: QuotientGroebnerResult

    @property
    def basis(self):
        """Return the certified generic lex Gröbner basis on this locus."""
        return self.certificate.basis

    def verify(self):
        """Check exact polynomial identities and the declared open guards."""
        if not self.certificate.verify():
            return False
        atoms = self.condition.args if self.condition.func is sp.And else (self.condition,)
        required = {
            *(sp.Ne(value, 0) for value in self.certificate.denominator_guards),
            *(sp.Eq(value, 0) for value in self.certificate.parameter_ideal.groebner_basis),
        }
        return required <= set(atoms)


@dataclass(frozen=True)
class PrimeParameterCover:
    """A partition into parameter loci with replayable generic bases."""

    pieces: tuple[PrimeCoverPiece, ...]
    parameters: tuple[sp.Symbol, ...]
    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]

    def verify(self):
        """Replay all proofs and verify symbolic coverage/disjointness."""
        return all(
            piece.verify()
            and piece.certificate.variables == self.variables
            and piece.certificate.equations == self.equations
            and piece.certificate.parameter_ideal.ring.variables == self.parameters
            for piece in self.pieces
        ) and check_parameter_partition(self.pieces, self.parameters)

    def specialize(self, values):
        """Evaluate the unique stratum basis at exact algebraic parameters."""
        values = tuple(sp.sympify(value) for value in values)
        if len(values) != len(self.parameters) or any(
            value.has(sp.Float) or value.is_algebraic is not True for value in values
        ):
            raise ValueError("one exact algebraic value is required for each parameter")
        substitutions = dict(zip(self.parameters, values, strict=True))
        matches = [
            piece
            for piece in self.pieces
            if sp.simplify(piece.condition.subs(substitutions, simultaneous=True)) is sp.true
        ]
        if len(matches) != 1:
            raise ArithmeticError("the supplied values do not select one certified stratum")
        return tuple(
            sp.cancel(expression.subs(substitutions, simultaneous=True))
            for expression in matches[0].basis
        )


def _nonvanishing_clauses(ideals):
    """Split each complement of an algebraic set into disjoint basic opens."""
    conditions = [()]
    for ideal in ideals:
        next_conditions = []
        polynomials = tuple(ideal.groebner_basis)
        for prefix in conditions:
            vanished = []
            for polynomial in polynomials:
                next_conditions.append((*prefix, *vanished, sp.Ne(polynomial, 0)))
                vanished.append(sp.Eq(polynomial, 0))
        conditions = next_conditions
    return conditions


def comprehensive_prime_cover(
    equations, variables, parameters, *, max_nodes=40, max_pieces=80, max_pairs=150
):
    """Recursively certify a Gröbner cover using prime quotient function fields.

    For each branch V(P) minus specified closed loci, compute a generic basis
    in Fr(QQ[a]/P). Every new inverted coefficient defines an exceptional
    variety. Recursion proceeds through its field-relative minimal primes;
    overlaps are partitioned by exact polynomial equations/inequations.
    Failure to certify an exceptional ideal raises instead of returning a
    partial cover. The result is bounded because recursive Gröbner covers can
    have exponential complexity.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    if not variables or not parameters or set(variables) & set(parameters):
        raise ValueError("disjoint nonempty fiber variables and parameters are required")
    if min(max_nodes, max_pieces, max_pairs) < 1:
        raise ValueError("resource limits must be positive")
    equations = tuple(sp.sympify(expression) for expression in equations)
    for expression in equations:
        sp.Poly(expression, *variables, *parameters, domain=sp.QQ)
    parameter_ring = PolynomialRing(parameters, sp.QQ, "lex")
    pieces = []
    nodes = 0

    def visit(prime, inequalities, extra):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise ParameterCoverLimitError("parameter branch budget exhausted")
        # On V(P), a previously inverted polynomial may have become zero.
        if any(prime.contains(expression) for expression in inequalities):
            return
        result = groebner_over_prime_parameters(equations, variables, prime, max_pairs=max_pairs)
        guards = tuple(
            guard
            for guard in result.denominator_guards
            if not any(sp.cancel(guard - previous) == 0 for previous in inequalities)
        )
        conditions = (
            *(sp.Eq(polynomial, 0) for polynomial in prime.groebner_basis),
            *(sp.Ne(value, 0) for value in (*inequalities, *guards)),
            *extra,
        )
        condition = sp.And(*conditions)
        if condition is not sp.false:
            if len(pieces) >= max_pieces:
                raise ParameterCoverLimitError("parameter-piece budget exhausted")
            pieces.append(PrimeCoverPiece(condition, result))
        # The locus of failure is the union of the hypersurfaces guard=0.
        # Earlier guards are excluded, making these branches disjoint.
        prior = list(inequalities)
        for guard in guards:
            exceptional = prime + parameter_ring.ideal(guard)
            if exceptional.is_unit:
                prior.append(guard)
                continue
            candidates = tuple(exceptional.minimal_primes())
            if not candidates:
                raise ParameterCoverLimitError("exceptional variety has no certified primes")
            for index, candidate in enumerate(candidates):
                # Disjointly partition V(candidate) minus earlier components.
                for clause in _nonvanishing_clauses(candidates[:index]):
                    visit(
                        candidate,
                        tuple(prior),
                        (*extra, *clause),
                    )
            prior.append(guard)

    visit(parameter_ring.ideal(), (), ())
    cover = PrimeParameterCover(tuple(pieces), parameters, variables, equations)
    if not cover.verify():
        raise ArithmeticError("recursive parameter cover failed symbolic verification")
    return cover


__all__ = [
    "ParameterCoverLimitError",
    "PrimeCoverPiece",
    "PrimeParameterCover",
    "comprehensive_prime_cover",
]
