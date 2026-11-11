"""Certified algebraic constants of finitely generated function fields.

For a prime ideal P over QQ, a rational function p/q in Frac(QQ[x]/P)
may satisfy a polynomial over QQ. Its minimal polynomial is obtained by
eliminating the coordinates from the saturated graph ideal
(P, t*q-p):q^infinity. A proposed constant field is accepted as complete
only when every conjugate graph component is certified geometrically prime
and their intersection equals the extended original ideal.

For a prime hypersurface quadratic in one coordinate, the square class of
its discriminant gives a complete, terminating constant-field algorithm:
nonconstant squarefree factors prove geometric integrality, while a constant
nonsquare square class gives a quadratic constant field. A standalone proof
replays the discriminant identity and, in the split case, the saturated graph
ideals. Other presentations still have a bounded candidate search; failure to
discover a suitable primitive constant never means the constant field is QQ.
"""

from dataclasses import dataclass
from itertools import product

import sympy as sp

from .ideals import Ideal, PolynomialRing


@dataclass(frozen=True)
class AlgebraicConstants:
    """A certified full constant field and its geometric generic components.

    The fraction ``numerator/denominator`` represents a primitive algebraic
    constant in the original function field. Its ``minimal_polynomial`` is
    irreducible over QQ. Every ``geometric_primes`` entry is proven absolutely
    prime over ``field``; exact intersection reconstructs the base-changed
    input ideal. When the constants are QQ, numerator is zero and denominator
    is one.
    """

    field: object
    splitting_field: object
    minimal_polynomial: sp.Expr
    numerator: sp.Expr
    denominator: sp.Expr
    geometric_primes: tuple[Ideal, ...]
    certificate: object = None

    def verify(self, prime):
        """Replay a complete discriminant certificate, when available.

        Discriminant-derived results carry an independently replayable
        certificate. Search results without that evidence return False.
        """
        return self.certificate is not None and self.certificate.verify(prime, self)


@dataclass(frozen=True)
class QuadraticConstantCertificate:
    """Recheckable completeness proof for a quadratic function field over QQ.

    If ``square_factor`` is None, the discriminant has a nonconstant factor
    with odd multiplicity, proving irreducibility in QQbar(t)[y]. Otherwise
    ``discriminant = square_class * square_factor**2``. Since the defining
    prime has degree two over QQ(t), its full constant field is QQ(sqrt(c)).
    All divisions and component graph constructions are independently
    reconstructed by the verifier; the producer's Gröbner basis is not trusted.
    """

    coordinate: sp.Symbol
    discriminant: sp.Expr
    square_class: sp.Rational | None
    square_factor: sp.Expr | None

    def verify(self, prime, result):
        """Independently recheck discriminant, full constants, and geometric graphs."""
        try:
            if not isinstance(prime, Ideal) or prime.ring.domain != sp.QQ:
                return False
            if len(prime.generators) != 1 or prime.is_unit:
                return False
            polynomial = sp.Poly(prime.generators[0], *prime.ring.variables, domain=sp.QQ)
            if sp.Poly(polynomial.as_expr(), self.coordinate).degree() != 2:
                return False
            factors = polynomial.factor_list()[1]
            if len(factors) != 1 or factors[0][1] != 1:
                return False
            parameters = tuple(v for v in prime.ring.variables if v != self.coordinate)
            quadratic = sp.Poly(
                polynomial.as_expr(), self.coordinate, domain=sp.QQ.frac_field(*parameters)
            )
            a, b, d = (sp.cancel(value) for value in quadratic.all_coeffs())
            delta = sp.cancel(b * b - 4 * a * d)
            if sp.cancel(delta - self.discriminant) != 0:
                return False
            if self.square_class is None:
                if self.square_factor is not None:
                    return False
                numerator, denominator = sp.fraction(delta)
                factors = sp.Poly(
                    sp.expand(numerator * denominator), *parameters, domain=sp.QQ
                ).factor_list()[1]
                if not any(
                    exponent % 2 for factor, exponent in factors if factor.total_degree() > 0
                ):
                    return False
                minimal = sp.Poly(result.minimal_polynomial)
                return (
                    minimal.degree() == 1
                    and minimal.nth(0) == 0
                    and minimal.LC() == 1
                    and result.field == sp.QQ
                    and result.splitting_field == sp.QQ
                    and result.geometric_primes == (prime,)
                    and result.numerator == 0
                    and result.denominator == 1
                )
            c = self.square_class
            if not isinstance(c, sp.Rational) or self.square_factor is None:
                return False
            if sp.cancel(delta - c * self.square_factor**2) != 0:
                return False
            if sp.sqrt(c).is_Rational:
                return False
            primitive = 2 * a * self.coordinate + b
            if sp.cancel(result.numerator - primitive) != 0:
                return False
            if sp.cancel(result.denominator - self.square_factor) != 0:
                return False
            if not prime.contains(result.numerator**2 - c * result.denominator**2):
                return False
            variable = next(iter(sp.Poly(result.minimal_polynomial).gens))
            if sp.Poly(result.minimal_polynomial, variable, domain=sp.QQ).monic() != sp.Poly(
                variable**2 - c, variable, domain=sp.QQ
            ):
                return False
            root = sp.sqrt(c)
            field = sp.QQ.algebraic_field(root)
            if result.field != field or result.splitting_field != field:
                return False
            if len(result.geometric_primes) != 2:
                return False
            from .ideals import PolynomialRing

            extended = PolynomialRing(prime.ring.variables, field, prime.ring.order)
            reference = prime.change_ring(extended)
            # These are graphs y=(-b +/- sqrt(c)*s)/(2a), localized at a.
            # Their contractions are geometrically prime by Gauss's lemma.
            actual = tuple(
                extended.ideal(
                    *prime.generators, primitive - sign * root * self.square_factor
                ).saturate(a)
                for sign in (1, -1)
            )
            if (
                actual[0] == actual[1]
                or any(
                    not any(candidate == expected for expected in actual)
                    for candidate in result.geometric_primes
                )
                or result.geometric_primes[0] == result.geometric_primes[1]
            ):
                return False
            if any(
                not _is_geometrically_prime_linear_graph(component, self.coordinate)
                for component in actual
            ):
                return False
            return actual[0].intersection(actual[1]) == reference
        except (
            ValueError,
            TypeError,
            ArithmeticError,
            sp.PolynomialError,
            NotImplementedError,
        ):
            return False


def _quadratic_discriminant_data(equation, coordinate, variables):
    """Return the exact square class of a quadratic over a rational-function field.

    Nonconstant odd squarefree factors stay non-squares after any algebraic
    constant extension in characteristic zero. Hence ``c=None`` is a complete
    geometric-integrality proof, not a bounded search result.
    """
    parameters = tuple(v for v in variables if v != coordinate)
    if not parameters:
        return None
    q = sp.Poly(equation, coordinate, domain=sp.QQ.frac_field(*parameters))
    if q.degree() != 2:
        return None
    a, b, d = (sp.cancel(v) for v in q.all_coeffs())
    delta = sp.cancel(b * b - 4 * a * d)
    numerator, denominator = sp.fraction(delta)
    # Multiplication by denominator**2 preserves the square class.
    polynomial = sp.Poly(sp.expand(numerator * denominator), *parameters, domain=sp.QQ)
    if polynomial.is_zero:
        return None  # characteristic-zero prime quadratics are separable
    content, factors = polynomial.sqf_list()
    if any(exponent % 2 for factor, exponent in factors if factor.total_degree() > 0):
        return a, b, d, delta, None, None
    root_part = sp.prod(factor.as_expr() ** (exponent // 2) for factor, exponent in factors)
    square_factor = sp.cancel(root_part / denominator)
    c = sp.Rational(content)
    if sp.cancel(delta - c * square_factor**2) != 0:
        raise ArithmeticError("discriminant squareclass reconstruction failed")
    return a, b, d, delta, c, square_factor


def _is_geometrically_prime_linear_graph(ideal, coordinate):
    """Prove a primitive linear equation in y remains prime over QQbar.

    For R=k[t], a*y+b is absolutely irreducible when gcd_R(a,b)=1:
    over kbar(t) it has degree one, and gcd remains one over kbar[t].
    """
    variables = ideal.ring.variables
    parameters = tuple(v for v in variables if v != coordinate)
    if not parameters or len(ideal.groebner_basis) != 1:
        return False
    polynomial = sp.Poly(
        ideal.groebner_basis[0], coordinate, domain=ideal.ring.domain.frac_field(*parameters)
    )
    if polynomial.degree() != 1:
        return False
    expression = ideal.groebner_basis[0]
    coeffs = sp.Poly(expression, coordinate).all_coeffs()
    if len(coeffs) != 2:
        return False
    first = sp.Poly(coeffs[0], *parameters, domain=ideal.ring.domain)
    second = sp.Poly(coeffs[1], *parameters, domain=ideal.ring.domain)
    return first.gcd(second).total_degree() == 0


def _quadratic_constants(prime):
    """Complete constant computation for a prime quadratic hypersurface."""
    if len(prime.generators) != 1 or prime.ring.domain != sp.QQ:
        return None
    variables = prime.ring.variables
    for coordinate in variables:
        data = _quadratic_discriminant_data(prime.generators[0], coordinate, variables)
        if data is None:
            continue
        a, b, _, delta, c, s = data
        certificate = QuadraticConstantCertificate(coordinate, delta, c, s)
        if c is None:
            symbol = sp.Dummy("constant")
            result = AlgebraicConstants(
                sp.QQ, sp.QQ, symbol, sp.S.Zero, sp.S.One, (prime,), certificate
            )
        else:
            if sp.sqrt(c).is_Rational:
                raise ValueError("a rational-square discriminant contradicts primality")
            from .ideals import PolynomialRing

            root = sp.sqrt(c)
            field = sp.QQ.algebraic_field(root)
            ambient = PolynomialRing(variables, field, prime.ring.order)
            primitive = sp.cancel(2 * a * coordinate + b)
            components = tuple(
                ambient.ideal(*prime.generators, primitive - sign * root * s).saturate(a)
                for sign in (1, -1)
            )
            symbol = sp.Dummy("constant")
            result = AlgebraicConstants(
                field, field, symbol**2 - c, primitive, s, components, certificate
            )
        if not result.verify(prime):
            raise ArithmeticError("quadratic constant-field proof failed independent replay")
        return result
    return None


def _prime_over_base(ideal):
    """Certify relative primality, without inferring it from reconstruction."""
    if ideal.is_unit:
        return False
    if ideal.is_zero:
        return True
    if len(ideal.generators) == 1:
        poly = sp.Poly(ideal.generators[0], *ideal.ring.variables, domain=ideal.ring.domain)
        factors = poly.factor_list()[1]
        return len(factors) == 1 and factors[0][1] == 1
    from .geometric import _graph_prime_decomposition

    if _graph_prime_decomposition(ideal) is not None:
        return True
    from .components import DecompositionLimitError, associated_primes

    try:
        return associated_primes(ideal) == (ideal,)
    except (DecompositionLimitError, NotImplementedError):
        return None


def _smooth_point(ideal, *, max_points=200):
    """Find a ground-field smooth point with a full-rank Jacobian witness.

    At a smooth rational point of an integral finite-type scheme, exactly one
    geometric component passes through the point. Galois transitivity of
    geometric components of a ground-field integral scheme proves geometric
    integrality. This is a sufficient test, not a smoothness requirement on
    the whole variety.
    """
    variables = ideal.ring.variables
    codimension = len(variables) - ideal.krull_dimension
    if codimension == 0:
        return (sp.S.Zero,) * len(variables)
    equations = ideal.generators
    jacobian = sp.Matrix([[sp.diff(f, v) for v in variables] for f in equations])
    candidates = [sp.S.Zero, sp.S.One, -sp.S.One, sp.S(2), sp.S(-2)]
    if getattr(ideal.ring.domain, "is_AlgebraicField", False):
        for generator in getattr(ideal.ring.domain, "orig_ext", ()):
            if generator not in candidates:
                candidates.extend((generator, -generator))
    for index, point in enumerate(product(candidates, repeat=len(variables))):
        if index >= max_points:
            break
        substitution = dict(zip(variables, point, strict=True))
        if any(sp.simplify(f.subs(substitution, simultaneous=True)) != 0 for f in equations):
            continue
        if jacobian.subs(substitution, simultaneous=True).rank() == codimension:
            return point
    return None


def certify_geometric_prime(ideal, *, max_points=200):
    """Prove a prime ideal stays prime over the algebraic closure.

    Supports arbitrary-dimensional nonlinear ideals when they have a smooth
    rational point and are provably prime over the coefficient field. Also
    recognizes linear and polynomial-graph primes. Raises NotImplementedError
    when the available criteria cannot certify geometric primality; it does
    not treat an inconclusive search as evidence of reducibility.
    """
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    if ideal.is_unit:
        return False
    if ideal.is_zero:
        return True
    from .geometric import _graph_prime_decomposition

    variables = ideal.ring.variables
    # A proper affine-linear ideal has a polynomial-ring quotient, so no
    # general primary-decomposition call is needed for its certificate.
    if all(sp.Poly(g, *variables).total_degree() <= 1 for g in ideal.groebner_basis):
        return True
    if _graph_prime_decomposition(ideal) is not None:
        return True
    if len(ideal.generators) == 1:
        for variable in variables:
            if _is_geometrically_prime_linear_graph(ideal, variable):
                return True
    if ideal.ring.domain == sp.QQ and len(ideal.generators) == 1 and len(variables) >= 2:
        from .geometric import _projectively_smooth

        polynomial = sp.Poly(ideal.generators[0], *variables, domain=sp.QQ)
        if polynomial.total_degree() <= 4 and _projectively_smooth(polynomial.as_expr(), variables):
            return True
    primality = _prime_over_base(ideal)
    if primality is None:
        raise NotImplementedError("field-relative primality could not be certified")
    if not primality:
        return False
    if ideal.ring.domain == sp.QQ and len(ideal.generators) == 1:
        for coordinate in variables:
            data = _quadratic_discriminant_data(ideal.generators[0], coordinate, variables)
            if data is not None:
                # A constant nonsquare is an obstruction to geometric
                # primality; an odd squarefree discriminant factor proves it.
                return data[4] is None
    if _smooth_point(ideal, max_points=max_points) is not None:
        return True
    raise NotImplementedError("no geometric-primality certificate found for the prime ideal")


def _constant_polynomial(prime, numerator, denominator, symbol, max_degree):
    variables = prime.ring.variables
    if prime.contains(denominator):
        return None
    extended = PolynomialRing((symbol, *variables), sp.QQ, "lex")
    graph = extended.ideal(*prime.generators, symbol * denominator - numerator)
    relation = graph.saturate(denominator).eliminate(*variables)
    if relation.is_zero or relation.is_unit:
        return None
    basis = relation.groebner_basis
    if len(basis) != 1:
        return None
    polynomial = sp.Poly(basis[0], symbol, domain=sp.QQ).monic()
    if not 1 < polynomial.degree() <= max_degree:
        return None
    factors = polynomial.factor_list()[1]
    if len(factors) != 1 or factors[0][1] != 1:
        raise ArithmeticError(
            "a constant in a prime function field has reducible minimal polynomial"
        )
    degree = polynomial.degree()
    identity = sum(
        polynomial.nth(power) * numerator**power * denominator ** (degree - power)
        for power in range(degree + 1)
    )
    if not prime.contains(identity):
        raise ArithmeticError("candidate constant failed its defining polynomial identity")
    return polynomial


def _rational_candidates(variables, degree, limit):
    pool = [*variables]
    if degree >= 2:
        pool.extend(x * y for index, x in enumerate(variables) for y in variables[index:])
    used = set()
    for numerator in pool:
        for denominator in (sp.S.One, *pool):
            if numerator == denominator:
                continue
            key = (numerator, denominator)
            if key in used:
                continue
            used.add(key)
            yield key
            if len(used) >= limit:
                return


def algebraic_constant_field(
    prime, *, max_degree=4, numerator_degree=1, max_candidates=60, max_points=200
):
    """Certify the full algebraic constant field of a supported QQ function field.

    Searches low-degree rational functions for a primitive algebraic constant.
    A candidate is accepted only when all conjugate graph ideals are certified
    geometrically prime and their intersection recovers the extended prime.
    Hence the result describes the *full* algebraic closure of QQ in the
    function field, not merely a known algebraic subfield.

    A geometrically prime field has constants QQ. Unsupported primes or a
    search that does not discover a full primitive constant raise
    ``NotImplementedError``. The bounded search is not a general decision
    procedure for all finitely generated function fields.
    """
    if not isinstance(prime, Ideal) or prime.ring.domain != sp.QQ:
        raise TypeError("a polynomial prime ideal over QQ is required")
    if prime.krull_dimension <= 0:
        raise ValueError("a positive-dimensional prime ideal is required")
    if min(max_degree, numerator_degree, max_candidates, max_points) < 1:
        raise ValueError("all search limits must be positive")
    primality = _prime_over_base(prime)
    if primality is None:
        raise NotImplementedError("field-relative primality could not be certified")
    if not primality:
        raise ValueError("the ideal is not certified prime over QQ")
    symbol = sp.Dummy("constant")
    quadratic = _quadratic_constants(prime)
    if quadratic is not None:
        return quadratic
    try:
        is_geometric = certify_geometric_prime(prime, max_points=max_points)
    except NotImplementedError:
        is_geometric = None
    if is_geometric is True:
        return AlgebraicConstants(sp.QQ, sp.QQ, symbol, sp.S.Zero, sp.S.One, (prime,))

    variables = prime.ring.variables
    for numerator, denominator in _rational_candidates(variables, numerator_degree, max_candidates):
        polynomial = _constant_polynomial(prime, numerator, denominator, symbol, max_degree)
        if polynomial is None:
            continue
        roots = sp.roots(polynomial.as_expr(), symbol)
        if len(roots) != polynomial.degree() or any(
            multiplicity != 1 for multiplicity in roots.values()
        ):
            continue
        extension = sp.QQ.algebraic_field(*roots)
        constant_field = sp.QQ.algebraic_field(next(iter(roots)))
        ambient = PolynomialRing(variables, extension, prime.ring.order)
        extended_prime = prime.change_ring(ambient)
        components = []
        for root in roots:
            component = ambient.ideal(*prime.generators, numerator - root * denominator)
            component = component.saturate(denominator)
            try:
                if not certify_geometric_prime(component, max_points=max_points):
                    break
            except NotImplementedError:
                break
            components.append(component)
        else:
            reconstructed = ambient.ideal(1)
            for component in components:
                reconstructed = reconstructed.intersection(component)
            if reconstructed != extended_prime:
                raise ArithmeticError("constant-field components failed reconstruction")
            return AlgebraicConstants(
                constant_field,
                extension,
                polynomial.as_expr(),
                numerator,
                denominator,
                tuple(components),
            )
    raise NotImplementedError("full algebraic constant field could not be certified")


__all__ = [
    "AlgebraicConstants",
    "QuadraticConstantCertificate",
    "algebraic_constant_field",
    "certify_geometric_prime",
]
