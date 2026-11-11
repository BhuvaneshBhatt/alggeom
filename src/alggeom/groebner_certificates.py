"""Exact certificates for polynomial bases supplied by algroots.

Basis construction and membership-witness production belong to algroots.
This module checks the resulting ideal identities, reducedness and S-pair
remainders before a certificate is used in geometric calculations.
"""

from collections.abc import Mapping
from dataclasses import dataclass

from algroots.domains import polynomial
from algroots.groebner_proofs import construct_membership_witnesses
from sympy.polys.groebnertools import groebner as reference_groebner
from sympy.polys.polyerrors import PolynomialError


def _divides(left, right):
    return all(a <= b for a, b in zip(left, right, strict=True))


def _poly(expression, ring):
    if getattr(expression, "ring", None) is ring:
        return expression
    return ring.from_dict(polynomial(expression, ring.symbols, ring.domain).rep.to_dict())


def _sparse_membership(candidate, row, generators, ring):
    """Replay a mixed sparse/expression row without expanding coefficient products."""
    accumulator = {}
    for factor, generator in zip(row, generators, strict=True):
        coefficients = (
            factor.items() if isinstance(factor, Mapping) else _poly(factor, ring).terms()
        )
        generator_terms = generator.terms()
        for powers, scalar in coefficients:
            powers = tuple(powers)
            if len(powers) != len(ring.symbols) or any(
                not isinstance(power, int) or power < 0 for power in powers
            ):
                return False
            scalar = ring.domain.convert(scalar)
            if not scalar:
                continue
            for exponents, generator_scalar in generator_terms:
                product = tuple(a + b for a, b in zip(powers, exponents, strict=True))
                value = accumulator.get(product, ring.domain.zero) + scalar * generator_scalar
                # Cancellation must remove support, so equality can compare exact sparse maps.
                if value:
                    accumulator[product] = value
                else:
                    accumulator.pop(product, None)
    return accumulator == dict(candidate.terms())


def _verify_memberships(generators, basis, witnesses, ring, qq_reducer):
    """Check witness shape and exact identities using the selected arithmetic."""
    rows = tuple(tuple(row) for row in witnesses)
    if len(rows) != len(basis) or any(len(row) != len(generators) for row in rows):
        return False
    all_sparse = all(all(isinstance(factor, Mapping) for factor in row) for row in rows)
    all_expanded = all(all(not isinstance(factor, Mapping) for factor in row) for row in rows)
    if qq_reducer is not None and (all_sparse or all_expanded):
        if all_sparse:
            return qq_reducer.verify_sparse_memberships(generators, basis, rows)
        converted = tuple(tuple(_poly(factor, ring) for factor in row) for row in rows)
        return qq_reducer.verify_expanded_memberships(generators, basis, converted)
    for candidate, row in zip(basis, rows, strict=True):
        if any(isinstance(factor, Mapping) for factor in row):
            if not _sparse_membership(candidate, row, generators, ring):
                return False
        else:
            reconstruction = ring.zero
            for factor, generator in zip(row, generators, strict=True):
                reconstruction += _poly(factor, ring) * generator
            if reconstruction != candidate:
                return False
    return True


def _verify_critical_pairs(basis, leading, ring, qq_reducer):
    """Apply Buchberger's criterion to a monic candidate with cached leading terms."""
    for index, left in enumerate(basis):
        for other_index, right in enumerate(basis[:index]):
            left_leading, right_leading = leading[index], leading[other_index]
            common = tuple(max(a, b) for a, b in zip(left_leading, right_leading, strict=True))
            # Relatively prime leading monomials satisfy the product criterion.
            if all(a + b == c for a, b, c in zip(left_leading, right_leading, common, strict=True)):
                continue
            if qq_reducer is not None:
                if qq_reducer.s_pair_remainder(index, other_index):
                    return False
            else:
                left_factor = ring.from_dict(
                    {tuple(m - a for m, a in zip(common, left_leading, strict=True)): 1}
                )
                right_factor = ring.from_dict(
                    {tuple(m - b for m, b in zip(common, right_leading, strict=True)): 1}
                )
                if (left_factor * left - right_factor * right).rem(basis):
                    return False
    return True


@dataclass(frozen=True)
class GroebnerBasisProof:
    """A reduced basis with exact generator-combination evidence.

    Each memberships row expresses one basis element using the original
    generators. Coefficients may be polynomial expressions or sparse rational
    maps. A compact_proof stores a rational arithmetic DAG instead of rows.
    Supply one representation and call verify to check its claims.
    """

    basis: tuple
    memberships: tuple
    compact_proof: object = None

    def verify(self, generators, ring, *, normal_form_backend="sympy"):
        """Replay the proof against the supplied original generators and ring."""
        if self.compact_proof is not None and self.memberships:
            return False
        return verify_reduced_groebner_basis(
            generators,
            self.basis,
            ring,
            membership_witnesses=(None if self.compact_proof is not None else self.memberships),
            compact_proof=self.compact_proof,
            normal_form_backend=normal_form_backend,
        )


def verify_reduced_groebner_basis(
    generators,
    candidate,
    ring,
    *,
    membership_witnesses=None,
    compact_proof=None,
    normal_form_backend="sympy",
):
    """Verify a reduced monic basis independently of its construction backend.

    Membership witnesses have one row per candidate and one column per input
    generator. Exact replay establishes candidate inclusion in the input
    ideal; without witnesses, an independent reference basis checks this
    inclusion. Input reduction establishes the reverse inclusion. Reducedness
    and every required S-pair remainder are checked in both paths.
    """
    if compact_proof is not None and membership_witnesses is not None:
        return False
    if normal_form_backend not in ("sympy", "flint"):
        raise ValueError("normal_form_backend must be sympy or flint")
    if normal_form_backend == "flint" and not getattr(ring.domain, "is_QQ", False):
        raise ValueError("flint rational normal form requires QQ")
    if not getattr(ring.domain, "is_Field", False):
        raise ValueError("Groebner certification requires a coefficient field")
    try:
        original = tuple(_poly(value, ring) for value in generators)
        basis = tuple(_poly(value, ring) for value in candidate)
        if any(not value or ring.domain.one != value.LC for value in basis):
            return False
        if not basis:
            if membership_witnesses is not None and tuple(membership_witnesses):
                return False
            return not any(original)
        leading = [tuple(value.LM) for value in basis]
        if len(set(leading)) != len(leading):
            return False
        if any(
            _divides(other, powers)
            for index, value in enumerate(basis)
            for powers, _ in value.terms()
            for other_index, other in enumerate(leading)
            if other_index != index
        ):
            return False
        qq_reducer = None
        if normal_form_backend == "flint":
            from algroots.qq_normal_form import FLINTQQNormalForm

            qq_reducer = FLINTQQNormalForm(basis, ring)
        if any(
            (qq_reducer.remainder(value) if qq_reducer is not None else value.rem(basis))
            for value in original
        ):
            return False
        if compact_proof is not None:
            if qq_reducer is not None:
                valid = qq_reducer.verify_rational_dag(compact_proof, original, basis)
            else:
                from algroots.compact_rational_proof import verify_rational_dag

                valid = verify_rational_dag(compact_proof, original, basis, ring)
            if not valid:
                return False
        elif membership_witnesses is None:
            independent_basis = tuple(reference_groebner(original, ring))
            if any(value.rem(independent_basis) for value in basis):
                return False
        elif not _verify_memberships(original, basis, membership_witnesses, ring, qq_reducer):
            return False
        return _verify_critical_pairs(basis, leading, ring, qq_reducer)
    except (ValueError, TypeError, ArithmeticError, AttributeError, PolynomialError):
        # Malformed coefficients and proof shapes are failed certificates.
        return False


__all__ = [
    "GroebnerBasisProof",
    "construct_membership_witnesses",
    "verify_reduced_groebner_basis",
]
