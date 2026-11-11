"""Exact parameter-partition proofs and replayable polynomial reduction identities.

The verifier uses ideal arithmetic on constructible loci, exact polynomial
identities, and Buchberger's criterion. It does not re-run the cover dispatcher.
"""

from dataclasses import dataclass

import sympy as sp

from .constructible import ConstructibleSet, LocallyClosedSet
from .ideals import PolynomialRing


class CoverProofError(ValueError):
    """A cover condition or certificate has no supported exact proof."""


def _conditions(condition):
    if condition is sp.true:
        return ()
    if condition is sp.false:
        return (sp.false,)
    if condition.func is sp.And:
        return condition.args
    return (condition,)


def _parameter_domain(parameters, conditions):
    if not parameters:
        return sp.QQ
    extensions = set()
    for condition in conditions:
        for atom in _conditions(condition):
            if atom in (sp.true, sp.false):
                continue
            if not isinstance(atom, (sp.Equality, sp.Unequality)):
                raise CoverProofError(f"unsupported parameter condition: {atom}")
            expression = sp.together(atom.lhs - atom.rhs)
            numerator, denominator = sp.fraction(expression)
            for value in (numerator, denominator):
                polynomial = sp.Poly(value, *parameters, extension=True)
                domain = polynomial.domain
                if getattr(domain, "is_AlgebraicField", False):
                    extensions.update(domain.orig_ext)
                elif domain in (sp.ZZ_I, sp.QQ_I):
                    extensions.add(sp.I)
                elif domain not in (sp.ZZ, sp.QQ):
                    raise CoverProofError("parameter conditions must have algebraic coefficients")
    return (
        sp.QQ.algebraic_field(*sorted(extensions, key=sp.default_sort_key)) if extensions else sp.QQ
    )


def _piece_as_constructible(condition, ring):
    equalities = []
    nonzero = []
    for atom in _conditions(condition):
        if atom is sp.false:
            return ConstructibleSet(())
        if atom is sp.true:
            continue
        if not isinstance(atom, (sp.Equality, sp.Unequality)):
            raise CoverProofError(f"unsupported parameter condition: {atom}")
        numerator, denominator = sp.fraction(sp.together(atom.lhs - atom.rhs))
        # Rational conditions are meaningful only where their denominators are defined.
        if denominator != 1:
            nonzero.append(denominator)
        if isinstance(atom, sp.Equality):
            equalities.append(numerator)
        else:
            nonzero.append(numerator)
    return ConstructibleSet((LocallyClosedSet(ring.ideal(*equalities), tuple(nonzero)),))


def _is_empty(locus):
    return all(piece.closure_ideal.is_unit for piece in locus.pieces)


def _propositional_partition(conditions, parameters):
    """Prove a partition from factor identities alone, when possible.

    Over a field, a product vanishes exactly when one factor vanishes. Replacing
    each irreducible polynomial factor by a Boolean variable preserves this
    implication. A propositional partition is therefore also a geometric
    partition. Algebraic dependencies between distinct factors may prevent
    this shortcut from proving a valid cover; the exact ideal test then runs.
    """
    if len(conditions) > 24 or len(parameters) > 4:
        return False
    atoms = {}

    def polynomial_condition(expression, zero):
        polynomial = sp.Poly(expression, *parameters, domain=sp.QQ_I)
        if polynomial.is_zero:
            return sp.true if zero else sp.false
        if polynomial.total_degree() == 0:
            return sp.false if zero else sp.true
        if polynomial.total_degree() > 8:
            raise ValueError("large factorization")
        factors = sp.factor_list(polynomial)[1]
        flags = []
        for factor, _ in factors:
            monic = factor.monic().as_expr()
            if monic not in atoms:
                atoms[monic] = sp.Dummy("factor_zero", boolean=True)
            flags.append(atoms[monic])
        return sp.Or(*flags) if zero else sp.And(*(sp.Not(flag) for flag in flags))

    def condition_formula(condition):
        if condition in (sp.true, sp.false):
            return condition
        pieces = []
        for atom in _conditions(condition):
            if atom in (sp.true, sp.false):
                pieces.append(atom)
                continue
            if not isinstance(atom, (sp.Equality, sp.Unequality)):
                raise ValueError("unsupported parameter condition")
            numerator, denominator = sp.fraction(sp.together(atom.lhs - atom.rhs))
            pieces.append(polynomial_condition(denominator, False))
            pieces.append(polynomial_condition(numerator, isinstance(atom, sp.Equality)))
        return sp.And(*pieces)

    try:
        formulas = tuple(condition_formula(c) for c in conditions)
        if len(atoms) > 12:
            return False
        if sp.simplify_logic(sp.Or(*formulas), force=True) is not sp.true:
            return False
        return all(
            sp.simplify_logic(sp.And(first, second), force=True) is sp.false
            for index, first in enumerate(formulas)
            for second in formulas[:index]
        )
    except (sp.PolynomialError, sp.CoercionFailed, ValueError, TypeError):
        return False


def check_parameter_partition(pieces, parameters):
    """Prove coverage and pairwise disjointness over an algebraic closure.

    The proof is set-theoretic and uses saturation, ideal intersections and
    exact constructible-set difference; it is not based on test points.
    """
    parameters = tuple(parameters)
    conditions = tuple(piece.condition for piece in pieces)
    domain = _parameter_domain(parameters, conditions)
    if not parameters:
        nonempty = sum(condition is not sp.false for condition in conditions)
        return nonempty == 1
    if _propositional_partition(conditions, parameters):
        return True
    ring = PolynomialRing(parameters, domain)
    loci = tuple(_piece_as_constructible(condition, ring) for condition in conditions)
    for i, first in enumerate(loci):
        for second in loci[:i]:
            if not _is_empty(first & second):
                return False
    remaining = ConstructibleSet((LocallyClosedSet(ring.ideal()),))
    for locus in loci:
        remaining = remaining - locus
    return _is_empty(remaining)


@dataclass(frozen=True)
class ReductionIdentity:
    """A polynomial equality dividend = sum(quotient*divisor) + remainder."""

    dividend: sp.Expr
    divisors: tuple[sp.Expr, ...]
    quotients: tuple[sp.Expr, ...]
    remainder: sp.Expr

    def verify(self):
        """Check the stored polynomial identity exactly over its coefficient field."""
        if len(self.divisors) != len(self.quotients):
            return False
        difference = (
            self.dividend
            - self.remainder
            - sum(
                (
                    quotient * divisor
                    for quotient, divisor in zip(self.quotients, self.divisors, strict=True)
                ),
                sp.S.Zero,
            )
        )
        return sp.cancel(difference) == 0


@dataclass(frozen=True)
class StratumProof:
    """An exact Gröbner-basis proof on a constructible parameter stratum."""

    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    basis: tuple[sp.Expr, ...]
    source: tuple[sp.Expr, ...]
    basis_membership: tuple[ReductionIdentity, ...]
    source_reductions: tuple[ReductionIdentity, ...]
    critical_pairs: tuple[ReductionIdentity, ...]
    nonzero_guards: tuple[sp.Expr, ...]
    parameter_bindings: tuple[tuple[sp.Symbol, sp.Expr], ...] = ()

    def verify(self):
        """Check all identities, leading terms and guards without Buchberger search."""
        if len(self.basis_membership) != len(self.basis):
            return False
        if len(self.source_reductions) != len(self.source):
            return False
        if not all(
            record.verify()
            for record in (*self.basis_membership, *self.source_reductions, *self.critical_pairs)
        ):
            return False
        if any(sp.cancel(record.remainder) != 0 for record in self.basis_membership):
            return False
        if any(
            sp.cancel(record.remainder) != 0 or record.divisors != self.basis
            for record in (*self.source_reductions, *self.critical_pairs)
        ):
            return False
        if any(
            record.divisors != self.source or record.dividend != basis
            for basis, record in zip(self.basis, self.basis_membership, strict=True)
        ):
            return False
        if len(self.critical_pairs) != len(self.basis) * (len(self.basis) - 1) // 2:
            return False
        index = 0
        for i, first in enumerate(self.basis):
            for second in self.basis[:i]:
                pair = self.critical_pairs[index]
                index += 1
                if pair.dividend != _s_polynomial(first, second, self.variables):
                    return False
        # A zero S-remainder certifies a Gröbner basis only if its quotient
        # terms obey the standard-representation leading-monomial bound.
        for expression in self.basis:
            if sp.cancel(sp.Poly(expression, *self.variables).LC(order="lex") - 1) != 0:
                return False
        position = 0
        for i, left in enumerate(self.basis):
            left_lm = sp.Poly(left, *self.variables).LM(order="lex").exponents
            for right in self.basis[:i]:
                right_lm = sp.Poly(right, *self.variables).LM(order="lex").exponents
                bound = tuple(max(a, b) for a, b in zip(left_lm, right_lm, strict=True))
                for quotient, divisor in zip(
                    self.critical_pairs[position].quotients, self.basis, strict=True
                ):
                    product = sp.Poly(quotient * divisor, *self.variables)
                    if not product.is_zero and product.LM(order="lex").exponents > bound:
                        return False
                position += 1
        required = _denominators(
            (
                *self.basis,
                *self.source,
                *(
                    quotient
                    for record in (
                        *self.basis_membership,
                        *self.source_reductions,
                        *self.critical_pairs,
                    )
                    for quotient in record.quotients
                ),
            ),
            self.variables,
            self.parameters,
        )
        if not required.issubset(set(self.nonzero_guards)):
            return False
        # The proof is performed after replacing bound parameters by algebraic
        # coordinates. Each replacement must hold everywhere on the claimed
        # stratum; otherwise the identities would certify a different fiber.
        try:
            domain = _parameter_domain(self.parameters, (self.condition,))
            ring = PolynomialRing(self.parameters, domain)
            locus = _piece_as_constructible(self.condition, ring)
            substitutions = {}
            for parameter, replacement in reversed(self.parameter_bindings):
                if parameter not in self.parameters or parameter in substitutions:
                    return False
                substitutions[parameter] = sp.cancel(replacement.subs(substitutions))
            for parameter, replacement in substitutions.items():
                discrepancy = sp.cancel(parameter - replacement)
                if not _is_empty(
                    locus & _piece_as_constructible(sp.Ne(discrepancy, 0, evaluate=False), ring)
                ):
                    return False
            for guard in self.nonzero_guards:
                if guard == 0:
                    return False
                if guard.free_symbols & set(self.parameters):
                    bad = locus & _piece_as_constructible(sp.Eq(guard, 0, evaluate=False), ring)
                    if not _is_empty(bad):
                        return False
        except (sp.PolynomialError, sp.CoercionFailed, ValueError, TypeError):
            return False
        return True


def _s_polynomial(first, second, variables):
    first_poly = sp.Poly(first, *variables)
    second_poly = sp.Poly(second, *variables)
    first_lm = first_poly.LM(order="lex").exponents
    second_lm = second_poly.LM(order="lex").exponents
    multiple = tuple(max(left, right) for left, right in zip(first_lm, second_lm, strict=True))
    first_factor = sp.prod(
        variable ** (power - exponent)
        for variable, power, exponent in zip(variables, multiple, first_lm, strict=True)
    )
    second_factor = sp.prod(
        variable ** (power - exponent)
        for variable, power, exponent in zip(variables, multiple, second_lm, strict=True)
    )
    return sp.expand(first_factor * first - second_factor * second)


def _denominators(expressions, variables, parameters):
    guards = set()
    for expression in expressions:
        polynomial = sp.Poly(expression, *variables) if variables else None
        coefficients = polynomial.coeffs() if polynomial is not None else (expression,)
        for coefficient in coefficients:
            _, denominator = sp.fraction(sp.cancel(coefficient))
            if denominator.free_symbols & set(parameters):
                guards.add(sp.factor(denominator))
    return guards


def prove_stratum(piece, equations, variables, parameters):
    """Construct identities using a traced Buchberger computation on one branch.

    The trace expresses every basis element in the original generators.
    Subsequent source and S-polynomial reductions certify ideal equality and
    Buchberger's criterion. All coefficient denominators are recorded.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    substitutions = {}
    for parameter, replacement in reversed(piece.parameter_bindings):
        substitutions[parameter] = sp.cancel(replacement.subs(substitutions))
    source = tuple(sp.cancel(eq.subs(substitutions, simultaneous=True)) for eq in equations)
    target = tuple(sp.cancel(eq.subs(substitutions, simultaneous=True)) for eq in piece.basis)
    free = tuple(piece.free_parameters)
    constants = tuple(piece.algebraic_constants)
    base = sp.QQ.algebraic_field(*constants) if constants else sp.QQ
    domain = base.frac_field(*free) if free else base
    count = len(source)
    basis = []
    representations = []
    guards = set()
    for index, expression in enumerate(source):
        poly = sp.Poly(expression, *variables, domain=domain)
        if poly.is_zero:
            continue
        leading = poly.LC(order="lex")
        guards.update(_denominators((1 / leading,), (), parameters))
        basis.append(poly.monic().as_expr())
        row = [sp.S.Zero] * count
        row[index] = 1 / leading
        representations.append(tuple(row))
    pairs = [(i, j) for i in range(len(basis)) for j in range(i)]
    while pairs:
        i, j = pairs.pop()
        first_lm = sp.Poly(basis[i], *variables, domain=domain).LM(order="lex").exponents
        second_lm = sp.Poly(basis[j], *variables, domain=domain).LM(order="lex").exponents
        multiple = tuple(max(a, b) for a, b in zip(first_lm, second_lm, strict=True))
        first_factor = sp.prod(
            v ** (m - a) for v, m, a in zip(variables, multiple, first_lm, strict=True)
        )
        second_factor = sp.prod(
            v ** (m - b) for v, m, b in zip(variables, multiple, second_lm, strict=True)
        )
        s_polynomial = first_factor * basis[i] - second_factor * basis[j]
        quotients, remainder = sp.reduced(
            s_polynomial, basis, *variables, order="lex", domain=domain
        )
        guards.update(_denominators((*quotients, remainder), variables, parameters))
        if sp.Poly(remainder, *variables, domain=domain).is_zero:
            continue
        row = [
            first_factor * representations[i][n]
            - second_factor * representations[j][n]
            - sum((quotients[k] * representations[k][n] for k in range(len(basis))), sp.S.Zero)
            for n in range(count)
        ]
        leading = sp.Poly(remainder, *variables, domain=domain).LC(order="lex")
        row = tuple(sp.cancel(value / leading) for value in row)
        guards.update(_denominators((*row, 1 / leading), variables, parameters))
        basis.append(sp.Poly(remainder, *variables, domain=domain).monic().as_expr())
        representations.append(row)
        pairs.extend((len(basis) - 1, k) for k in range(len(basis) - 1))
    # The independent trace need not use the cover's exact basis order, but it
    # must reproduce the same specialized ideal and leading monomial ideal.
    actual = sp.groebner(target, *variables, order="lex", domain=domain)
    generated = sp.groebner(basis, *variables, order="lex", domain=domain)
    if tuple(actual.polys) != tuple(generated.polys):
        raise CoverProofError("the recorded branch basis disagrees with the Buchberger trace")
    # If the trace and cover bases differ, trace-basis membership witnesses
    # cannot be reused directly for cover elements; refuse an unsupported proof.
    if tuple(basis) != target:
        raise CoverProofError("the branch basis order differs from its traced construction")
    basis_membership = tuple(
        ReductionIdentity(value, source, row, sp.S.Zero)
        for value, row in zip(basis, representations, strict=True)
    )
    source_reductions = []
    for expression in source:
        q, r = (
            sp.reduced(expression, basis, *variables, order="lex", domain=domain)
            if basis
            else ([], expression)
        )
        if sp.cancel(r) != 0:
            raise CoverProofError("an original generator is not generated by the recorded basis")
        source_reductions.append(
            ReductionIdentity(
                expression, target, tuple(q) + (sp.S.Zero,) * (len(target) - len(q)), sp.S.Zero
            )
        )
        guards.update(_denominators(q, variables, parameters))
    critical_pairs = []
    for i, first in enumerate(target):
        for second in target[:i]:
            s = _s_polynomial(first, second, variables)
            q, r = sp.reduced(s, target, *variables, order="lex", domain=domain)
            if sp.cancel(r) != 0:
                raise CoverProofError("an S-polynomial does not reduce to zero")
            critical_pairs.append(
                ReductionIdentity(
                    s, target, tuple(q) + (sp.S.Zero,) * (len(target) - len(q)), sp.S.Zero
                )
            )
            guards.update(_denominators(q, variables, parameters))
    guards.update(_denominators((*target, *source), variables, parameters))
    proof = StratumProof(
        piece.condition,
        variables,
        parameters,
        target,
        source,
        basis_membership,
        tuple(source_reductions),
        tuple(critical_pairs),
        tuple(sorted(guards, key=sp.default_sort_key)),
        tuple(piece.parameter_bindings),
    )
    if not proof.verify():
        raise CoverProofError("a reduction identity or coefficient guard failed verification")
    return proof


@dataclass(frozen=True)
class GroebnerCoverProof:
    """A parameter partition with exact polynomial proofs for every basis."""

    pieces: tuple[StratumProof, ...]
    parameters: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]

    def verify(self):
        """Replay source identities, reductions and the exact parameter partition."""
        for proof in self.pieces:
            if proof.variables != self.variables or proof.parameters != self.parameters:
                return False
            substitutions = {}
            for parameter, replacement in reversed(proof.parameter_bindings):
                substitutions[parameter] = sp.cancel(replacement.subs(substitutions))
            expected = tuple(
                sp.cancel(eq.subs(substitutions, simultaneous=True)) for eq in self.equations
            )
            if len(proof.source) != len(expected) or any(
                sp.cancel(found - original) != 0
                for found, original in zip(proof.source, expected, strict=True)
            ):
                return False
            if not proof.verify():
                return False
        return check_parameter_partition(self.pieces, self.parameters)


__all__ = [
    "CoverProofError",
    "ReductionIdentity",
    "StratumProof",
    "GroebnerCoverProof",
    "check_parameter_partition",
    "prove_stratum",
]
