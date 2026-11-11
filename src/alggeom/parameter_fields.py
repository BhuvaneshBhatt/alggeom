"""Exact fraction fields of prime polynomial parameter quotients.

A prime ideal P in QQ[a] defines the integral domain A=QQ[a]/P and its
fraction field Frac(A). Field elements are fractions of polynomial normal
forms, with equality tested by cross multiplication modulo P. This makes
rational arithmetic possible without solving P for a particular parameter.
"""

from dataclasses import dataclass

import sympy as sp

from .ideals import Ideal


class ParameterFieldError(ValueError):
    """The supplied parameter ideal has no certified fraction field."""


def _prime_certificate(ideal):
    """Establish primeness before permitting fraction-field inversion.

    The fast paths are mathematically sufficient: a nonconstant irreducible
    polynomial generates a prime ideal over a field, and polynomial-graph
    quotient rings are polynomial rings. The fallback uses associated primes
    and checks that the ideal *itself* is prime, not merely primary.
    """
    if ideal.is_zero:
        return True
    if ideal.is_unit:
        return False
    if len(ideal.generators) == 1:
        expression = ideal.generators[0]
        polynomial = sp.Poly(expression, *ideal.ring.variables, domain=sp.QQ)
        factors = polynomial.factor_list()[1]
        return len(factors) == 1 and factors[0][1] == 1
    from .geometric import _graph_prime_decomposition

    if _graph_prime_decomposition(ideal) is not None:
        return True
    from .components import DecompositionLimitError, associated_primes

    try:
        return associated_primes(ideal) == (ideal,)
    except (DecompositionLimitError, NotImplementedError):
        return False


class PrimeParameterField:
    """Fraction field of ``QQ[parameters]/prime_ideal``.

    The constructor proves that the ideal is prime before allowing division;
    a primary but nonprime ideal is rejected. Computations use exact Gröbner
    normal forms and work on implicit nonlinear parameter varieties without
    selecting a radical or rational parametrization.
    """

    def __init__(self, prime_ideal: Ideal):
        if not isinstance(prime_ideal, Ideal) or prime_ideal.ring.domain != sp.QQ:
            raise TypeError("a polynomial ideal over QQ is required")
        if not prime_ideal.ring.variables:
            raise ValueError("at least one parameter is required")
        if not _prime_certificate(prime_ideal):
            raise ParameterFieldError("the parameter ideal is not certified prime")
        self.prime_ideal = prime_ideal
        self.parameters = prime_ideal.ring.variables
        self.inverted = set()

    def reduce(self, expression):
        """Return the canonical polynomial normal form modulo the prime ideal."""
        expression = sp.sympify(expression)
        sp.Poly(expression, *self.parameters, domain=sp.QQ)
        return sp.expand(self.prime_ideal.normal_form(expression))

    def element(self, numerator, denominator=1):
        """Construct a quotient-field element, rejecting zero divisors."""
        numerator = self.reduce(numerator)
        denominator = self.reduce(denominator)
        if denominator == 0:
            raise ZeroDivisionError("denominator vanishes in the parameter quotient")
        return ParameterFraction(self, numerator, denominator)

    def coerce(self, value):
        """Convert a rational expression or compatible field element."""
        if isinstance(value, ParameterFraction):
            if value.field is not self:
                raise TypeError("cannot mix distinct parameter fraction fields")
            return value
        value = sp.sympify(value)
        top, bottom = sp.fraction(sp.cancel(value))
        return self.element(top, bottom)

    def records(self):
        """List polynomial conditions required by prior field inversions."""
        return tuple(sorted(self.inverted, key=sp.default_sort_key))


@dataclass(frozen=True, eq=False)
class ParameterFraction:
    """An exact fraction in a certified prime parameter quotient field."""

    field: PrimeParameterField
    numerator: sp.Expr
    denominator: sp.Expr

    def _combine(self, other, operation):
        other = self.field.coerce(other)
        a, b = self.numerator, self.denominator
        c, d = other.numerator, other.denominator
        if operation == "add":
            return self.field.element(a * d + b * c, b * d)
        if operation == "multiply":
            return self.field.element(a * c, b * d)
        raise ValueError("unknown quotient-field operation")

    def __add__(self, other):
        return self._combine(other, "add")

    __radd__ = __add__

    def __neg__(self):
        return self.field.element(-self.numerator, self.denominator)

    def __sub__(self, other):
        return self + -self.field.coerce(other)

    def __rsub__(self, other):
        return self.field.coerce(other) - self

    def __mul__(self, other):
        return self._combine(other, "multiply")

    __rmul__ = __mul__

    def inverse(self):
        """Invert a nonzero residue class and record its specialization guard."""
        if self.numerator == 0:
            raise ZeroDivisionError("zero has no inverse")
        self.field.inverted.add(self.numerator)
        return self.field.element(self.denominator, self.numerator)

    def __truediv__(self, other):
        return self * self.field.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.field.coerce(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise TypeError("fraction-field powers must be integral")
        if exponent < 0:
            return self.inverse() ** -exponent
        result = self.field.element(1)
        base, power = self, exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power //= 2
        return result

    def __eq__(self, other):
        other = self.field.coerce(other)
        return (
            self.field.reduce(
                self.numerator * other.denominator - other.numerator * self.denominator
            )
            == 0
        )

    def as_expr(self):
        """Return a representative rational expression in the parameters."""
        return sp.cancel(self.numerator / self.denominator)


def _lead(poly):
    return max(poly) if poly else None


def _monomial_divides(divisor, dividend):
    return all(a <= b for a, b in zip(divisor, dividend, strict=True))


def _polynomial_add(first, second):
    result = first.copy()
    for monomial, coefficient in second.items():
        combined = result.get(monomial, coefficient.field.element(0)) + coefficient
        if combined == 0:
            result.pop(monomial, None)
        else:
            result[monomial] = combined
    return result


def _polynomial_shift(poly, monomial, factor):
    if factor == 0:
        return {}
    return {
        tuple(a + b for a, b in zip(power, monomial, strict=True)): factor * coefficient
        for power, coefficient in poly.items()
    }


def _polynomial_product(first, second, field):
    result = {}
    for first_power, first_coefficient in first.items():
        for second_power, second_coefficient in second.items():
            power = tuple(a + b for a, b in zip(first_power, second_power, strict=True))
            result = _polynomial_add(result, {power: first_coefficient * second_coefficient})
    return result


def _divide(poly, divisors, field):
    """Multivariate lex division over the quotient function field."""
    remainder, current = {}, poly.copy()
    quotients = [{} for _ in divisors]
    while current:
        head = _lead(current)
        leading = current[head]
        for index, divisor in enumerate(divisors):
            dhead = _lead(divisor)
            if dhead is not None and _monomial_divides(dhead, head):
                shift = tuple(a - b for a, b in zip(head, dhead, strict=True))
                scalar = leading / divisor[dhead]
                quotients[index] = _polynomial_add(quotients[index], {shift: scalar})
                negated = _polynomial_shift(divisor, shift, -scalar)
                current = _polynomial_add(current, negated)
                break
        else:
            remainder[head] = leading
            current = _polynomial_add(current, {head: -leading})
    return quotients, remainder


def _from_expression(expression, variables, field):
    poly = sp.Poly(expression, *variables)
    return {
        exponent: field.coerce(coefficient)
        for exponent, coefficient in poly.terms()
        if field.coerce(coefficient) != 0
    }


def _to_expression(poly, variables):
    return sp.Add(
        *(
            coefficient.as_expr()
            * sp.prod(v**power for v, power in zip(variables, powers, strict=True))
            for powers, coefficient in sorted(poly.items(), reverse=True)
        )
    )


@dataclass(frozen=True)
class QuotientGroebnerResult:
    """A generic lex basis with independently replayable membership witnesses.

    Identities hold in the fraction field of the declared prime quotient.
    Every coefficient inverse used by the construction is recorded as a
    nonvanishing guard. Guard conditions must be imposed on any parameter
    stratum before specializing the generic basis.
    """

    basis: tuple[sp.Expr, ...]
    denominator_guards: tuple[sp.Expr, ...]
    parameter_ideal: Ideal
    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    membership_witnesses: tuple[tuple[sp.Expr, ...], ...]

    def verify(self):
        """Verify source membership, Buchberger's criterion and all guards.

        The verifier never calls the basis constructor. Polynomial identities
        are checked modulo the parameter prime using exact fraction-field
        arithmetic; the quotient field must itself be certified prime.
        """
        field = PrimeParameterField(self.parameter_ideal)
        try:
            original = tuple(
                _from_expression(expression, self.variables, field) for expression in self.equations
            )
            basis = tuple(
                _from_expression(expression, self.variables, field) for expression in self.basis
            )
            if len(basis) != len(self.membership_witnesses):
                return False
            for polynomial, row in zip(basis, self.membership_witnesses, strict=True):
                if len(row) != len(original):
                    return False
                candidate = {}
                for coefficient, generator in zip(row, original, strict=True):
                    candidate = _polynomial_add(
                        candidate,
                        _polynomial_product(
                            _from_expression(coefficient, self.variables, field),
                            generator,
                            field,
                        ),
                    )
                if _polynomial_add(candidate, {key: -value for key, value in polynomial.items()}):
                    return False
            if any(_divide(generator, basis, field)[1] for generator in original):
                return False
            for index, first in enumerate(basis):
                for second in basis[:index]:
                    left, right = _lead(first), _lead(second)
                    if left is None or right is None:
                        return False
                    multiple = tuple(max(a, b) for a, b in zip(left, right, strict=True))
                    first_shift = tuple(a - b for a, b in zip(multiple, left, strict=True))
                    second_shift = tuple(a - b for a, b in zip(multiple, right, strict=True))
                    s_poly = _polynomial_add(
                        _polynomial_shift(first, first_shift, field.element(1) / first[left]),
                        _polynomial_shift(second, second_shift, -field.element(1) / second[right]),
                    )
                    if _divide(s_poly, basis, field)[1]:
                        return False
            # A missing guard is a material error: it can make a purported
            # generic identity undefined at an included parameter value.
            declared = {field.reduce(value) for value in self.denominator_guards}
            necessary = {field.reduce(value) for value in field.records()}
            for expression in (
                *self.basis,
                *(value for row in self.membership_witnesses for value in row),
            ):
                for coefficient in _from_expression(expression, self.variables, field).values():
                    necessary.add(field.reduce(coefficient.denominator))
            necessary.discard(sp.S.One)
            return sp.S.Zero not in declared and necessary <= declared
        except (ZeroDivisionError, ValueError, TypeError, sp.PolynomialError):
            return False


def groebner_over_prime_parameters(equations, variables, parameter_ideal, *, max_pairs=200):
    """Compute a lex Gröbner basis over ``Frac(QQ[a]/P)``.

    ``P`` must be a certified prime ideal. Every inverted residue numerator
    and coefficient denominator is returned as a nonvanishing guard; the
    result is a generic fiber basis on V(P) away from those guards. It makes
    no claim about fibers on the vanishing locus of a guard.
    """
    field = PrimeParameterField(parameter_ideal)
    variables = tuple(variables)
    if not variables or set(variables) & set(field.parameters):
        raise ValueError("distinct nonempty fiber variables are required")
    if max_pairs < 1:
        raise ValueError("max_pairs must be positive")
    equations = tuple(sp.sympify(eq) for eq in equations)
    for equation in equations:
        sp.Poly(equation, *variables, *field.parameters, domain=sp.QQ)
    original = [_from_expression(eq, variables, field) for eq in equations]
    basis = []
    witnesses = []
    zero_monomial = (0,) * len(variables)
    for index, poly in enumerate(original):
        if poly:
            lead = _lead(poly)
            scalar = poly[lead].inverse()
            basis.append(_polynomial_shift(poly, zero_monomial, scalar))
            row = [{} for _ in original]
            row[index] = {zero_monomial: scalar}
            witnesses.append(row)
    pairs = [(i, j) for i in range(len(basis)) for j in range(i)]
    count = 0
    while pairs:
        count += 1
        if count > max_pairs:
            raise RuntimeError("quotient-field Buchberger pair budget exceeded")
        i, j = pairs.pop()
        first, second = basis[i], basis[j]
        left, right = _lead(first), _lead(second)
        multiple = tuple(max(a, b) for a, b in zip(left, right, strict=True))
        first_shift = tuple(a - b for a, b in zip(multiple, left, strict=True))
        second_shift = tuple(a - b for a, b in zip(multiple, right, strict=True))
        s_poly = _polynomial_add(
            _polynomial_shift(first, first_shift, field.element(1)),
            _polynomial_shift(second, second_shift, field.element(-1)),
        )
        quotients, reduced = _divide(s_poly, basis, field)
        if reduced:
            row = []
            for source_index in range(len(original)):
                combination = _polynomial_add(
                    _polynomial_shift(witnesses[i][source_index], first_shift, field.element(1)),
                    _polynomial_shift(witnesses[j][source_index], second_shift, field.element(-1)),
                )
                for position, quotient in enumerate(quotients):
                    combination = _polynomial_add(
                        combination,
                        {
                            k: -v
                            for k, v in _polynomial_product(
                                quotient, witnesses[position][source_index], field
                            ).items()
                        },
                    )
                row.append(combination)
            scalar = reduced[_lead(reduced)].inverse()
            basis.append(_polynomial_shift(reduced, zero_monomial, scalar))
            witnesses.append([_polynomial_shift(poly, zero_monomial, scalar) for poly in row])
            pairs.extend((len(basis) - 1, k) for k in range(len(basis) - 1))
    # Every rational coefficient needs to remain defined under specialization.
    guards = set(field.records())
    for poly in (*basis, *(part for row in witnesses for part in row)):
        for coefficient in poly.values():
            guards.add(coefficient.denominator)
    guards = {field.reduce(guard) for guard in guards} - {sp.S.One}
    result = QuotientGroebnerResult(
        tuple(_to_expression(poly, variables) for poly in basis),
        tuple(sorted(guards, key=sp.default_sort_key)),
        parameter_ideal,
        variables,
        equations,
        tuple(tuple(_to_expression(poly, variables) for poly in row) for row in witnesses),
    )
    if not result.verify():
        raise ArithmeticError("quotient-field Gröbner certificate failed replay")
    return result


__all__ = [
    "ParameterFieldError",
    "PrimeParameterField",
    "ParameterFraction",
    "QuotientGroebnerResult",
    "groebner_over_prime_parameters",
]
