"""Automatic algebraic decomposition and candidate coverage accounting."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from algroots.curves import CurveError
from algroots.errors import QuotientAlgebraError
from sympy.polys.polyerrors import BasePolynomialError

from ._validation import (
    DecompositionError,
    constant_field,
    parse_relation,
    validate_integer,
    validate_symbols,
)
from .algebraic import AlgebraicReconstruction, normalize_factor, reconstruct_algebraic
from .factors import FactorDiscovery, discover_factors


@dataclass(frozen=True)
class CandidateAttempt:
    """A reconstruction result or an explicitly unresolved candidate diagnostic."""

    index: int
    result: AlgebraicReconstruction | None
    error: str | None = None


@dataclass(frozen=True)
class AlgebraicSearch:
    """Automatic minimal coefficient-field search, with exact coverage accounting."""

    discovery: FactorDiscovery
    parameter: sp.Symbol
    attempts: tuple[CandidateAttempt, ...]

    @property
    def relation(self):
        """The normalized input relation."""
        return self.discovery.relation

    @property
    def candidate_count(self):
        """Number of factor orbits discovered for reconstruction."""
        return len(self.discovery.orbits)

    @property
    def examined(self):
        """Number of candidate reconstruction attempts retained."""
        return len(self.attempts)

    @property
    def decompositions(self):
        """Every successful composition among the retained attempts."""
        return tuple(
            attempt.result.decomposition
            for attempt in self.attempts
            if attempt.result is not None and attempt.result.decomposition is not None
        )

    @property
    def complete(self):
        """Whether every relevant factor degree and orbit has a resolved reconstruction."""
        degree = self.relation.degree(self.relation.gens[0])
        return (
            self.discovery.complete
            and self.discovery.requested_degrees == tuple(range(2, degree + 1))
            and (self.examined == self.candidate_count)
            and all(attempt.result is not None for attempt in self.attempts)
        )

    @property
    def status(self):
        """Decomposability within the declared minimal coefficient-field scope."""
        if self.decompositions:
            return "decomposable"
        if not self.complete:
            return "incomplete"
        if self.relation.degree(self.relation.gens[0]) < 2:
            return "degree_below_two"
        return "indecomposable"

    def verify(self):
        """Replay orbit coverage and every claimed answer/rejection.

        Unresolved diagnostic text is not a proof of unsupportedness; it always
        prevents completeness and an indecomposability claim.
        """
        try:
            if not self.discovery.verify() or self.examined > self.candidate_count:
                return False
            validate_symbols((*self.relation.gens, self.parameter))
            if tuple(attempt.index for attempt in self.attempts) != tuple(range(self.examined)):
                return False
            for attempt, orbit in zip(self.attempts, self.discovery.orbits, strict=False):
                if attempt.result is None:
                    if not isinstance(attempt.error, str) or not attempt.error:
                        return False
                    continue
                reconstruction = attempt.result
                if (
                    attempt.error is not None
                    or reconstruction.relation != self.relation
                    or reconstruction.ambient is not orbit.ambient
                    or (
                        reconstruction.candidate
                        != normalize_factor(orbit.factor, self.relation.gens[0], orbit.ambient)
                    )
                    or reconstruction.parameter != self.parameter
                    or not reconstruction.verify()
                ):
                    return False
            return True
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


def search_algebraic(
    relation,
    variables,
    *,
    parameter=None,
    constants=sp.QQ,
    factor_degrees=None,
    max_degree=6,
    max_factor_dimension=32,
    max_separating_attempts=64,
    max_candidates=64,
    exhaustive=False,
    max_steps=256,
    max_dimension=256,
):
    """Search factor orbits for minimal coefficient-field compositions, retaining coverage evidence."""
    parameter = sp.Dummy("intermediate") if parameter is None else parameter
    constants = constant_field(constants)
    validate_integer(max_candidates, "max_candidates", 0)
    validate_integer(max_steps, "max_steps", 1)
    validate_integer(max_dimension, "max_dimension", 1)
    if not isinstance(exhaustive, bool):
        raise DecompositionError("exhaustive must be a Boolean")
    x, z = tuple(variables)
    validate_symbols((x, z, parameter))
    F = parse_relation(relation, (x, z), constants)
    degrees = tuple(range(2, F.degree(x) + 1)) if factor_degrees is None else tuple(factor_degrees)
    if any(d == 1 for d in degrees):
        raise DecompositionError("degree-one factors only give trivial minimal outer components")
    discovery = discover_factors(
        F,
        (x, z),
        degrees=degrees,
        constants=constants,
        max_degree=max_degree,
        max_factor_dimension=max_factor_dimension,
        max_separating_attempts=max_separating_attempts,
    )
    attempts = []
    for index, orbit in enumerate(discovery.orbits):
        if index >= max_candidates:
            break
        try:
            result = reconstruct_algebraic(
                F,
                (x, z),
                candidate=orbit.factor,
                ambient=orbit.ambient,
                parameter=parameter,
                constants=constants,
                max_steps=max_steps,
                max_dimension=max_dimension,
            )
        except (CurveError, QuotientAlgebraError, NotImplementedError) as exc:
            attempts.append(CandidateAttempt(index, None, f"{type(exc).__name__}: {exc}"))
            continue
        attempts.append(CandidateAttempt(index, result))
        if result.decomposition is not None and (not exhaustive):
            break
    return AlgebraicSearch(discovery, parameter, tuple(attempts))


__all__ = ["CandidateAttempt", "AlgebraicSearch", "search_algebraic"]
