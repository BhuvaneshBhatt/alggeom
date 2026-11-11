"""Absolute factorization of bivariate polynomials using logarithmic derivatives.

For a squarefree polynomial f in QQ[x, y], the space of closed rational
one-forms with poles along f and bounded bidegrees has one dimension per
geometric irreducible factor. Its multiplication algebra separates factors
using exact rational linear algebra and algebraic eigenvalues. Every factor
is recovered by polynomial gcd and checked by exact reconstruction.
"""

from dataclasses import dataclass

import sympy as sp


class FactorizationLimitError(NotImplementedError):
    """Exact factor reconstruction exceeds a declared algebraic complexity limit."""


@dataclass(frozen=True)
class AbsoluteFactors:
    """Geometrically irreducible factors and their multiplicities over a number field."""

    factors: tuple[tuple[sp.Expr, int], ...]
    field: object


def _closed_forms(polynomial, variables, max_unknowns):
    """Solve the logarithmic-derivative equations with bounded numerators."""
    x, y = variables
    degree_x = sp.degree(polynomial, x)
    degree_y = sp.degree(polynomial, y)
    if min(degree_x, degree_y) < 1:
        raise ValueError("the polynomial must involve both variables")
    exponents_x = [(i, j) for i in range(degree_x) for j in range(degree_y + 1)]
    exponents_y = [(i, j) for i in range(degree_x + 1) for j in range(degree_y)]
    total = len(exponents_x) + len(exponents_y)
    if total > max_unknowns:
        raise FactorizationLimitError("logarithmic form system exceeds max_unknowns")
    coefficients = sp.symbols(f"form_0:{total}")
    numerator_x = sum(
        coefficient * x**i * y**j
        for coefficient, (i, j) in zip(coefficients, exponents_x, strict=False)
    )
    numerator_y = sum(
        coefficient * x**i * y**j
        for coefficient, (i, j) in zip(coefficients[len(exponents_x) :], exponents_y, strict=False)
    )
    closed = sp.Poly(
        sp.expand(
            polynomial * (sp.diff(numerator_x, y) - sp.diff(numerator_y, x))
            - numerator_x * sp.diff(polynomial, y)
            + numerator_y * sp.diff(polynomial, x)
        ),
        x,
        y,
    )
    matrix, _ = sp.linear_eq_to_matrix(closed.coeffs(), coefficients)
    return tuple(
        (
            sp.expand(numerator_x.subs(dict(zip(coefficients, vector, strict=True)))),
            sp.expand(numerator_y.subs(dict(zip(coefficients, vector, strict=True)))),
        )
        for vector in matrix.nullspace()
    )


def _multiplication_matrices(polynomial, variables, forms):
    """Construct the étale residue algebra of the logarithmic-derivative space."""
    x, y = variables
    for shift in range(len(forms) + 2):
        derivative = sp.diff(polynomial, x) + shift * sp.diff(polynomial, y)
        if sp.gcd(polynomial, derivative) == 1:
            break
    else:
        raise ArithmeticError("no derivative separates all geometric branches")
    numerators = tuple(first + shift * second for first, second in forms)
    groebner = sp.groebner((polynomial,), x, y, order="lex", domain=sp.QQ)
    remainder_cache = {}

    def reduced(expression):
        terms = sp.Poly(expression, x, y).terms()
        result = 0
        for exponents, coefficient in terms:
            if exponents not in remainder_cache:
                monomial = x ** exponents[0] * y ** exponents[1]
                remainder_cache[exponents] = groebner.reduce(monomial)[1]
            result += coefficient * remainder_cache[exponents]
        return sp.Poly(result, x, y)

    dimension = len(forms)
    unknowns = sp.symbols(f"residue_0:{dimension}")
    candidate = derivative * sum(
        unknown * numerator for unknown, numerator in zip(unknowns, numerators, strict=True)
    )
    matrices = []
    for left in numerators:
        columns = []
        for right in numerators:
            equations = reduced(left * right - candidate).coeffs()
            solutions = sp.linsolve(equations, unknowns)
            if solutions is sp.EmptySet:
                raise ArithmeticError("logarithmic residue multiplication is not closed")
            answer = next(iter(solutions))
            if any(value.free_symbols & set(unknowns) for value in answer):
                raise ArithmeticError("logarithmic residue multiplication is not unique")
            columns.append(answer)
        matrices.append(sp.Matrix.hstack(*(sp.Matrix(column) for column in columns)))
    return tuple(matrices), numerators, derivative


def _squarefree_factors(polynomial, variables, max_unknowns, max_roots):
    forms = _closed_forms(polynomial, variables, max_unknowns)
    dimension = len(forms)
    if dimension == 1:
        return (polynomial,), sp.QQ
    if dimension < 1 or dimension > max_roots:
        raise FactorizationLimitError("absolute factor count exceeds max_roots")
    matrices, numerators, derivative = _multiplication_matrices(polynomial, variables, forms)
    candidates = (
        tuple(range(1, dimension + 1)),
        tuple(2**i for i in range(dimension)),
        tuple(i**2 + 1 for i in range(dimension)),
    )
    for weights in candidates:
        operator = sp.zeros(dimension)
        for weight, matrix in zip(weights, matrices, strict=True):
            operator += weight * matrix
        characteristic = sp.Poly(operator.charpoly().as_expr(), domain=sp.QQ)
        if sp.gcd(characteristic, characteristic.diff()).degree() == 0:
            break
    else:
        raise FactorizationLimitError("no separating element found in the residue algebra")
    roots = sp.roots(characteristic.as_expr(), characteristic.gen)
    if len(roots) != dimension or any(multiplicity != 1 for multiplicity in roots.values()):
        raise FactorizationLimitError("algebraic residue eigenvalues cannot be represented exactly")
    field = (
        sp.QQ.algebraic_field(*roots)
        if any(value.is_Rational is not True for value in roots)
        else sp.QQ
    )
    numerator = sum(weight * value for weight, value in zip(weights, numerators, strict=True))
    x, y = variables
    exact = sp.Poly(polynomial, x, y, domain=field)
    factors = []
    for eigenvalue in roots:
        polynomial_value = sp.Poly(numerator - eigenvalue * derivative, x, y, domain=field)
        factor = sp.gcd(exact, polynomial_value).monic()
        if factor.total_degree() < 1:
            raise ArithmeticError("a residue eigenvalue did not recover an algebraic component")
        factors.append(factor.as_expr())
    reconstructed = sp.Poly(sp.prod(factors), x, y, domain=field)
    quotient = exact.exquo(reconstructed)
    if quotient.total_degree() != 0 or quotient.is_zero:
        raise ArithmeticError("absolute factors did not reconstruct the original polynomial")
    return tuple(factors), field


def absolute_factor_bivariate(expression, variables, *, max_unknowns=100, max_roots=12):
    """Factor a bivariate QQ polynomial into certified geometric factors.

    Multiplicities are preserved. Factors are represented over a common exact
    number field. The method uses logarithmic derivatives and characteristic-
    zero residue-algebra splitting; it does not rely on factorization over a
    generic rational-function field, which could split a geometrically prime
    polynomial incorrectly.
    """
    variables = tuple(variables)
    if len(variables) != 2 or any(not isinstance(var, sp.Symbol) for var in variables):
        raise ValueError("two polynomial variables are required")
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    if polynomial.total_degree() <= 0:
        raise ValueError("a nonconstant polynomial is required")
    factors = []
    extensions = []
    for rational_factor, multiplicity in polynomial.factor_list()[1]:
        factor = rational_factor.as_expr()
        if any(sp.degree(factor, coordinate) == 0 for coordinate in variables):
            # A univariate polynomial may require algebraic roots beyond radicals.
            coordinate = next(var for var in variables if sp.degree(factor, var) > 0)
            roots = sp.roots(factor, coordinate)
            if sum(roots.values()) != sp.degree(factor, coordinate):
                raise FactorizationLimitError("a univariate factor has unresolved algebraic roots")
            factors.extend(
                (coordinate - root, multiplicity * power) for root, power in roots.items()
            )
            extensions.extend(root for root in roots if root.is_Rational is not True)
            continue
        geometric_factors, field = _squarefree_factors(factor, variables, max_unknowns, max_roots)
        factors.extend((part, multiplicity) for part in geometric_factors)
        extensions.extend(getattr(field, "orig_ext", ()))
    domain = sp.QQ.algebraic_field(*extensions) if extensions else sp.QQ
    product = sp.Poly(
        sp.prod(part**multiplicity for part, multiplicity in factors),
        *variables,
        domain=domain,
    )
    residual = sp.Poly(expression, *variables, domain=domain).exquo(product)
    if residual.total_degree() != 0 or residual.is_zero:
        raise ArithmeticError("the geometric factorization failed reconstruction")
    return AbsoluteFactors(tuple(factors), domain)


def absolute_factor_hypersurface(expression, variables, *, max_unknowns=90, max_roots=12):
    """Factor a rational-coefficient multivariate hypersurface geometrically.

    For two variables the optimized bivariate backend is used. In higher
    dimensions, solve the closed-form equations for logarithmic differentials
    in all coordinates, then reconstruct geometric factors from residue
    eigenvalues and exact multivariate gcds. Only characteristic zero and
    rational coefficients are supported; expensive systems raise explicitly.
    """
    from itertools import product

    variables = tuple(variables)
    if len(variables) == 2:
        return absolute_factor_bivariate(
            expression, variables, max_unknowns=max_unknowns, max_roots=max_roots
        )
    if len(variables) < 2:
        raise ValueError("at least two variables are required")
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    if polynomial.total_degree() <= 0:
        raise ValueError("a nonconstant polynomial is required")

    def split_squarefree(factor):
        degrees = tuple(sp.degree(factor, variable) for variable in variables)
        supports = []
        for coordinate in range(len(variables)):
            if not degrees[coordinate]:
                supports.append(())
                continue
            ranges = (
                range(degree + 1 - int(index == coordinate)) for index, degree in enumerate(degrees)
            )
            supports.append(tuple(product(*ranges)))
        unknown_count = sum(len(support) for support in supports)
        if unknown_count > max_unknowns:
            raise FactorizationLimitError("multivariate residue system exceeds max_unknowns")
        unknowns = sp.symbols(f"form_0:{unknown_count}")
        numerators = []
        offset = 0
        for support in supports:
            numerator = sum(
                unknowns[offset + position]
                * sp.prod(
                    variable**power for variable, power in zip(variables, exponents, strict=True)
                )
                for position, exponents in enumerate(support)
            )
            numerators.append(sp.sympify(numerator))
            offset += len(support)
        equations = []
        for first in range(len(variables)):
            for second in range(first + 1, len(variables)):
                closed = sp.expand(
                    factor
                    * (
                        sp.diff(numerators[first], variables[second])
                        - sp.diff(numerators[second], variables[first])
                    )
                    - numerators[first] * sp.diff(factor, variables[second])
                    + numerators[second] * sp.diff(factor, variables[first])
                )
                equations.extend(sp.Poly(closed, *variables).coeffs())
        matrix, _ = sp.linear_eq_to_matrix(equations, unknowns)
        forms = tuple(
            tuple(
                sp.expand(poly.subs(dict(zip(unknowns, vector, strict=True))))
                for poly in numerators
            )
            for vector in matrix.nullspace()
        )
        count = len(forms)
        if count == 1:
            return (factor,), sp.QQ
        if count < 1 or count > max_roots:
            raise FactorizationLimitError("absolute factor count exceeds max_roots")
        for shift in range(1, count + 4):
            derivative = sum(
                shift**index * sp.diff(factor, variable) for index, variable in enumerate(variables)
            )
            if sp.gcd(factor, derivative) == 1:
                break
        else:
            raise ArithmeticError("no derivative is nonzero on every factor")
        values = tuple(
            sum(shift**index * component for index, component in enumerate(form)) for form in forms
        )
        basis = sp.groebner((factor,), *variables, order="lex", domain=sp.QQ)
        remainder_cache = {}

        def reduced(expression):
            result = 0
            for exponents, coefficient in sp.Poly(expression, *variables).terms():
                if exponents not in remainder_cache:
                    monomial = sp.prod(
                        symbol**power for symbol, power in zip(variables, exponents, strict=True)
                    )
                    remainder_cache[exponents] = basis.reduce(monomial)[1]
                result += coefficient * remainder_cache[exponents]
            return sp.Poly(result, *variables)

        coefficients = sp.symbols(f"residue_0:{count}")
        candidate = derivative * sum(
            coefficient * value for coefficient, value in zip(coefficients, values, strict=True)
        )
        matrices = []
        for left in values:
            columns = []
            for right in values:
                solutions = sp.linsolve(reduced(left * right - candidate).coeffs(), coefficients)
                if solutions is sp.EmptySet:
                    raise ArithmeticError("residue algebra multiplication failed")
                answer = next(iter(solutions))
                if any(value.free_symbols & set(coefficients) for value in answer):
                    raise ArithmeticError("residue algebra product is not unique")
                columns.append(answer)
            matrices.append(sp.Matrix.hstack(*(sp.Matrix(column) for column in columns)))
        for weights in (
            tuple(range(1, count + 1)),
            tuple(2**index for index in range(count)),
            tuple(index**2 + 1 for index in range(count)),
        ):
            operator = sum(
                (weight * matrix for weight, matrix in zip(weights, matrices, strict=True)),
                sp.zeros(count),
            )
            characteristic = sp.Poly(operator.charpoly().as_expr(), domain=sp.QQ)
            if sp.gcd(characteristic, characteristic.diff()).degree() == 0:
                break
        else:
            raise FactorizationLimitError("residue algebra has no separating element")
        roots = sp.roots(characteristic.as_expr(), characteristic.gen)
        if len(roots) != count or any(multiplicity != 1 for multiplicity in roots.values()):
            raise FactorizationLimitError("residue eigenvalues cannot be represented exactly")
        field = (
            sp.QQ.algebraic_field(*roots)
            if any(root.is_Rational is not True for root in roots)
            else sp.QQ
        )
        numerator = sum(weight * value for weight, value in zip(weights, values, strict=True))
        exact = sp.Poly(factor, *variables, domain=field)
        factors = tuple(
            sp.gcd(
                exact,
                sp.Poly(numerator - root * derivative, *variables, domain=field),
            )
            .monic()
            .as_expr()
            for root in roots
        )
        reconstructed = sp.Poly(sp.prod(factors), *variables, domain=field)
        if exact.exquo(reconstructed).total_degree() != 0:
            raise ArithmeticError("geometric hypersurface factors failed reconstruction")
        if any(sp.Poly(part, *variables).total_degree() <= 0 for part in factors):
            raise ArithmeticError("residue eigenvalue failed to identify a geometric factor")
        return factors, field

    result = []
    extensions = []
    for factor, power in polynomial.factor_list()[1]:
        geometric, field = split_squarefree(factor.as_expr())
        result.extend((part, power) for part in geometric)
        extensions.extend(getattr(field, "orig_ext", ()))
    field = sp.QQ.algebraic_field(*extensions) if extensions else sp.QQ
    combined = sp.Poly(sp.prod(part**power for part, power in result), *variables, domain=field)
    if sp.Poly(expression, *variables, domain=field).exquo(combined).total_degree() != 0:
        raise ArithmeticError("absolute hypersurface reconstruction failed")
    return AbsoluteFactors(tuple(result), field)


__all__ = [
    "AbsoluteFactors",
    "FactorizationLimitError",
    "absolute_factor_bivariate",
    "absolute_factor_hypersurface",
]
