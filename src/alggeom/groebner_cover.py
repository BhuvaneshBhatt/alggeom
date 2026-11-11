"""Certified parameter-specialized Groebner systems over exact coefficient fields.

The computation branches on leading coefficients encountered by Buchberger's
algorithm. An open branch retains a basis whose reductions remain valid after
specialization. Exceptional fibers are recomputed, not evaluated from a
potentially invalid generic basis. One-parameter families admit a complete
finite decomposition when the exceptional algebraic roots can be represented.
For multiple parameters, exceptional divisors with polynomial algebraic
parameter-graph charts can be recursively split. More general implicit
parameter varieties raise an explicit coverage error.
"""

from dataclasses import dataclass

import sympy as sp


class GroebnerCoverLimitError(NotImplementedError):
    """The parameter branching requires an unsupported algebraic specialization."""


@dataclass(frozen=True)
class GroebnerCoverPiece:
    """A constructible parameter condition and a valid Groebner basis on it."""

    condition: sp.Expr
    basis: tuple[sp.Expr, ...]
    order: str
    parameter_bindings: tuple[tuple[sp.Symbol, sp.Expr], ...] = ()
    free_parameters: tuple[sp.Symbol, ...] = ()
    algebraic_constants: tuple[sp.Expr, ...] = ()


@dataclass(frozen=True)
class GroebnerCover:
    """Certified basis selection for parameter-specialized polynomial ideals."""

    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    pieces: tuple[GroebnerCoverPiece, ...]

    def verify_partition(self):
        """Prove exact coverage and pairwise disjointness of parameter strata."""
        from .cover_proofs import check_parameter_partition

        return check_parameter_partition(self.pieces, self.parameters)

    def certificate(self):
        """Produce replayable ideal-membership and Buchberger certificates.

        Each proof stores polynomial reduction quotients, S-pair reductions,
        and coefficient-denominator guards checked on its entire stratum.
        """
        from .cover_proofs import GroebnerCoverProof, prove_stratum

        result = GroebnerCoverProof(
            tuple(
                prove_stratum(piece, self.equations, self.variables, self.parameters)
                for piece in self.pieces
            ),
            self.parameters,
            self.equations,
            self.variables,
        )
        if not result.verify():
            raise ArithmeticError("Gröbner-cover proof does not replay")
        return result

    def specialize(self, values):
        """Return a specialized basis at an exact algebraic parameter tuple."""
        values = tuple(sp.sympify(value) for value in values)
        if len(values) != len(self.parameters) or any(
            value.is_algebraic is not True for value in values
        ):
            raise ValueError("one exact algebraic value is required per parameter")
        substitutions = dict(zip(self.parameters, values, strict=True))
        matched = []
        for piece in self.pieces:
            evaluated = piece.condition.subs(substitutions, simultaneous=True)
            if evaluated is sp.true or sp.simplify(evaluated) is sp.true:
                matched.append(piece)
        if len(matched) != 1:
            raise ArithmeticError("the cover does not select exactly one basis")
        return tuple(
            sp.cancel(polynomial.subs(substitutions, simultaneous=True))
            for polynomial in matched[0].basis
        )

    def verify_at(self, values):
        """Independently compare the selected basis with a specialized Groebner basis."""
        values = tuple(sp.sympify(value) for value in values)
        substitutions = dict(zip(self.parameters, values, strict=True))
        candidate = self.specialize(values)
        originals = tuple(eq.subs(substitutions, simultaneous=True) for eq in self.equations)
        domain = (
            sp.QQ.algebraic_field(*values)
            if any(value.is_Rational is not True for value in values)
            else sp.QQ
        )
        expected = sp.groebner(originals, *self.variables, order="lex", domain=domain)
        actual = sp.groebner(candidate, *self.variables, order="lex", domain=domain)
        return tuple(expected.polys) == tuple(actual.polys)


def _branch_basis(
    equations, variables, parameters, algebraic_constants, max_pairs, *, domain_override=None
):
    """Run Buchberger over the generic coefficient field and collect guards.

    Every division uses a leading coefficient that must remain nonzero. By
    preserving these leading coefficients, plus all rational denominators,
    specialization commutes with each step of the recorded reduction history.
    The critical-pair criterion then certifies the resulting basis on the
    corresponding open parameter set.
    """
    base = sp.QQ.algebraic_field(*algebraic_constants) if algebraic_constants else sp.QQ
    domain = (
        domain_override
        if domain_override is not None
        else (base.frac_field(*parameters) if parameters else base)
    )
    basis = []
    guard_factors = []

    def guard(polynomial):
        for term in (
            polynomial.all_coeffs()
            if len(variables) == 1
            else (coefficient for _, coefficient in polynomial.terms())
        ):
            _, denominator = sp.cancel(term).as_numer_denom()
            if denominator.free_symbols & set(parameters):
                guard_factors.append(denominator)
        leading = polynomial.LC(order="lex")
        numerator = sp.cancel(leading).as_numer_denom()[0]
        if numerator.free_symbols & set(parameters):
            guard_factors.append(numerator)

    for equation in equations:
        polynomial = sp.Poly(equation, *variables, domain=domain)
        if polynomial.is_zero:
            continue
        guard(polynomial)
        basis.append(polynomial.monic())
    pairs = [(i, j) for i in range(len(basis)) for j in range(i)]
    processed = 0
    while pairs:
        processed += 1
        if processed > max_pairs:
            raise GroebnerCoverLimitError("Buchberger critical-pair budget exhausted")
        i, j = pairs.pop()
        first, second = basis[i], basis[j]
        first_lm = first.LM(order="lex").exponents
        second_lm = second.LM(order="lex").exponents
        multiple = tuple(max(a, b) for a, b in zip(first_lm, second_lm, strict=True))
        first_factor = sp.prod(
            variable ** (power - exponent)
            for variable, power, exponent in zip(variables, multiple, first_lm, strict=True)
        )
        second_factor = sp.prod(
            variable ** (power - exponent)
            for variable, power, exponent in zip(variables, multiple, second_lm, strict=True)
        )
        s_polynomial = first_factor * first.as_expr() - second_factor * second.as_expr()
        _, remainder = sp.reduced(
            s_polynomial, [p.as_expr() for p in basis], *variables, order="lex", domain=domain
        )
        result = sp.Poly(remainder, *variables, domain=domain)
        if result.is_zero:
            continue
        guard(result)
        basis.append(result.monic())
        pairs.extend((len(basis) - 1, index) for index in range(len(basis) - 1))
    guarded = sp.Integer(1)
    for expression in guard_factors:
        _, denominator = sp.cancel(expression).as_numer_denom()
        numerator = sp.cancel(expression).as_numer_denom()[0]
        for part in (numerator, denominator):
            if part.free_symbols & set(parameters):
                guarded *= part
    if parameters and domain_override is None:
        guarded = sp.Poly(guarded, *parameters, domain=base).sqf_part().as_expr()
    return tuple(p.as_expr() for p in basis), guarded


def _polynomial_exceptional_roots(polynomial, selected, remaining, max_degree):
    """Find all roots polynomial in the remaining parameters over a number field.

    A root involving radicals of free parameters does not define a globally
    regular chart. Such exceptional varieties need quotient-function-field
    computations, so this helper declines them rather than treating a chosen
    radical branch as a complete parameter cover.
    """
    degree = sp.degree(polynomial, selected)
    if degree <= 1 or degree > max_degree:
        return None
    roots = sp.roots(polynomial, selected)
    if sum(roots.values()) != degree:
        return None
    constants = set()
    for root in roots:
        try:
            coefficient_domain = sp.Poly(root, *remaining, extension=True).domain
        except (sp.PolynomialError, sp.CoercionFailed, ValueError, TypeError):
            return None
        if coefficient_domain in (sp.ZZ_I, sp.QQ_I):
            constants.add(sp.I)
        elif getattr(coefficient_domain, "is_AlgebraicField", False):
            constants.update(coefficient_domain.orig_ext)
        elif coefficient_domain not in (sp.QQ, sp.ZZ):
            return None
        if sp.cancel(polynomial.subs(selected, root)) != 0:
            return None
    return tuple(roots), tuple(sorted(constants, key=sp.default_sort_key))


def comprehensive_groebner_cover(
    equations,
    variables,
    parameters,
    *,
    max_pairs=150,
    max_pieces=64,
    max_exceptional_degree=6,
):
    """Build specialization-valid polynomial bases for tractable parameter strata.

    A complete answer is returned for univariate parameter spaces when all
    exceptional roots admit exact representations, and for multivariate
    parameter spaces when all exceptional branches admit polynomial parameter charts over an algebraic
    constant extension, or a supported generic finite extension. Unsupported
    implicit parameter varieties raise instead of returning an incomplete cover.
    """
    variables, parameters = tuple(variables), tuple(parameters)
    if not variables or len(set((*variables, *parameters))) != len(variables) + len(parameters):
        raise ValueError("distinct fiber and parameter symbols are required")
    equations = tuple(sp.sympify(expression) for expression in equations)
    if any(
        expression.has(sp.Float) or expression.free_symbols - set((*variables, *parameters))
        for expression in equations
    ):
        raise ValueError("exact polynomial equations in the declared variables are required")
    for expression in equations:
        sp.Poly(expression, *variables, *parameters, domain=sp.QQ)
    if max_pairs < 1 or max_pieces < 1 or max_exceptional_degree < 1:
        raise ValueError("limits must be positive")

    pieces = []

    def add_piece(piece):
        if len(pieces) >= max_pieces:
            raise GroebnerCoverLimitError("parameter-piece budget exhausted")
        pieces.append(piece)

    def search(current_equations, free_parameters, substitutions, conditions, algebraic_constants):
        if len(pieces) >= max_pieces:
            raise GroebnerCoverLimitError("parameter-piece budget exhausted")
        basis, guard = _branch_basis(
            current_equations, variables, free_parameters, algebraic_constants, max_pairs
        )
        mapped = tuple(sp.cancel(polynomial.xreplace(substitutions)) for polynomial in basis)
        base_condition = sp.And(*conditions)
        if guard == 0:
            raise ArithmeticError("generic leading coefficient product vanished")
        if not free_parameters or not guard.free_symbols:
            if base_condition is not sp.false:
                add_piece(
                    GroebnerCoverPiece(
                        base_condition,
                        mapped,
                        "lex",
                        tuple(substitutions.items()),
                        free_parameters,
                        algebraic_constants,
                    )
                )
            return
        # The discriminant of this Buchberger computation defines exactly the
        # exceptional divisor where one of its leading coefficients may drop.
        generic_condition = sp.And(base_condition, sp.Ne(guard.xreplace(substitutions), 0))
        if generic_condition is not sp.false:
            add_piece(
                GroebnerCoverPiece(
                    generic_condition,
                    mapped,
                    "lex",
                    tuple(substitutions.items()),
                    free_parameters,
                    algebraic_constants,
                )
            )
        earlier_nonzero = []
        for exceptional, _ in sp.factor_list(guard, *free_parameters)[1]:
            polynomial = exceptional.as_expr()
            candidates = [
                parameter for parameter in free_parameters if polynomial.free_symbols <= {parameter}
            ]
            if not candidates:
                # A linear relation permits an exact rational parametrization.
                candidates = [
                    parameter
                    for parameter in free_parameters
                    if sp.degree(polynomial, parameter) == 1
                ]
            if not candidates:
                # A nonlinear exceptional divisor may split into regular
                # algebraic graphs after extending the constant field.
                # Recurse on each graph, including intersections of graphs;
                # later graphs exclude previous ones to keep the cover disjoint.
                for selected in free_parameters:
                    remaining = tuple(p for p in free_parameters if p != selected)
                    if not remaining:
                        continue
                    resolved = _polynomial_exceptional_roots(
                        polynomial, selected, remaining, max_exceptional_degree
                    )
                    if resolved is None:
                        continue
                    roots, constants = resolved
                    previous_graphs = []
                    for root in roots:
                        reduced_equations = tuple(
                            sp.cancel(eq.subs(selected, root)) for eq in current_equations
                        )
                        substitutions_next = {
                            **substitutions,
                            selected: root.xreplace(substitutions),
                        }
                        equation = sp.Eq(selected, root.xreplace(substitutions))
                        search(
                            reduced_equations,
                            remaining,
                            substitutions_next,
                            (*conditions, *earlier_nonzero, *previous_graphs, equation),
                            (*algebraic_constants, *constants),
                        )
                        previous_graphs.append(sp.Ne(selected, root.xreplace(substitutions)))
                    earlier_nonzero.append(sp.Ne(polynomial.xreplace(substitutions), 0))
                    break
                else:
                    # Over the generic point of an irreducible parameter
                    # hypersurface, adjoining one coordinate yields a finite
                    # algebraic extension of a rational-function field. When
                    # Buchberger uses no additional parameter inverses, this
                    # basis specializes along the entire hypersurface.
                    from sympy.polys.agca.extensions import FiniteExtension

                    selected = next(
                        (
                            parameter
                            for parameter in free_parameters
                            if sp.degree(polynomial, parameter) > 0
                        ),
                        None,
                    )
                    remaining = tuple(param for param in free_parameters if param != selected)
                    if selected is None or not remaining:
                        raise GroebnerCoverLimitError(
                            "exceptional parameter hypersurface cannot be represented"
                        )
                    base_field = sp.QQ.frac_field(*remaining)
                    defining = sp.Poly(polynomial, selected, domain=base_field)
                    irreducible_factors = defining.factor_list()[1]
                    if len(irreducible_factors) != 1 or irreducible_factors[0][1] != 1:
                        raise GroebnerCoverLimitError(
                            "exceptional parameter variety is reducible over the generic field"
                        )
                    extension = FiniteExtension(defining.monic())
                    exceptional_basis, exceptional_guard = _branch_basis(
                        current_equations,
                        variables,
                        free_parameters,
                        (),
                        max_pairs,
                        domain_override=extension,
                    )
                    if exceptional_guard.free_symbols & set(free_parameters):
                        raise GroebnerCoverLimitError(
                            f"a secondary exceptional divisor on {polynomial} requires "
                            "nested parameter-ideal specialization"
                        )
                    condition = sp.And(
                        *conditions,
                        *earlier_nonzero,
                        sp.Eq(polynomial.xreplace(substitutions), 0),
                    )
                    if condition is not sp.false:
                        add_piece(
                            GroebnerCoverPiece(
                                condition,
                                tuple(poly.xreplace(substitutions) for poly in exceptional_basis),
                                "lex",
                                tuple(substitutions.items()),
                                free_parameters,
                                algebraic_constants,
                            )
                        )
                earlier_nonzero.append(sp.Ne(polynomial.xreplace(substitutions), 0))
                continue
            selected = candidates[0]
            degree = sp.degree(polynomial, selected)
            if degree > max_exceptional_degree:
                raise GroebnerCoverLimitError("exceptional root degree exceeds limit")
            if polynomial.free_symbols <= {selected}:
                roots = sp.roots(polynomial, selected)
                if sum(roots.values()) != degree:
                    raise GroebnerCoverLimitError("exceptional roots are not represented exactly")
            else:
                roots = {sp.cancel(sp.solve(polynomial, selected)[0]): 1}
                coefficient = sp.diff(polynomial, selected)
                if coefficient.free_symbols:
                    raise GroebnerCoverLimitError(
                        "linear exceptional coefficient vanishes on another parameter stratum"
                    )
            for root in roots:
                reduced_parameters = tuple(p for p in free_parameters if p != selected)
                reduced_equations = tuple(
                    sp.cancel(eq.subs(selected, root)) for eq in current_equations
                )
                substitution = {**substitutions, selected: root.xreplace(substitutions)}
                root_condition = sp.Eq(selected, root.xreplace(substitutions))
                if sp.And(*conditions, *earlier_nonzero, root_condition) is sp.false:
                    continue
                search(
                    reduced_equations,
                    reduced_parameters,
                    substitution,
                    (*conditions, *earlier_nonzero, root_condition),
                    (*algebraic_constants, root)
                    if root.is_algebraic is True and root.is_Rational is not True
                    else algebraic_constants,
                )
            earlier_nonzero.append(sp.Ne(polynomial.xreplace(substitutions), 0))

    search(equations, parameters, {}, (), ())
    return GroebnerCover(equations, variables, parameters, tuple(pieces))


__all__ = [
    "GroebnerCover",
    "GroebnerCoverPiece",
    "GroebnerCoverLimitError",
    "comprehensive_groebner_cover",
]
