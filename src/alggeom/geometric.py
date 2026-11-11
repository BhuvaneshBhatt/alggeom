"""Geometric integrality and exact primary factors over algebraic extensions.

Absolute algorithms handle hypersurfaces, bounded finite algebras, and
complete intersections in disjoint variable blocks over QQ. Unrestricted
positive-dimensional nonprincipal splitting is not supplied.
"""

from dataclasses import dataclass

import sympy as sp

from ._algebra import integer
from .ideals import Ideal


@dataclass(frozen=True)
class GeometricComponents:
    """Absolute components explicitly represented over a splitting number field."""

    components: tuple[Ideal, ...]
    field: object


def _quadric_rank(polynomial, variables):
    homogenizer = sp.Dummy("homogenizer")
    homogeneous = sp.Poly(polynomial, *variables).homogenize(homogenizer).as_expr()
    coordinates = (*variables, homogenizer)
    matrix = sp.hessian(homogeneous, coordinates) / 2
    return int(matrix.rank())


def is_geometrically_integral(ideal):
    """Decide geometric integrality for principal linear/quadratic hypersurfaces.

    For an affine plane conic, homogenizing the equation gives a ternary
    quadratic form. Over an algebraic closure its zero set is irreducible
    exactly when that form has rank three. A univariate polynomial describes
    one geometrically integral point precisely when its degree is one.

    Supported finite zero-dimensional ideals are also tested by splitting
    their geometric primary components and checking the reduced quotient.
    """
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    if ideal.is_unit:
        return False
    if ideal.is_zero:
        return True
    if len(ideal.generators) > 1 and ideal.ring.domain == sp.QQ:
        separated = _separated_hypersurface_primary_decomposition(ideal)
        if separated is not None:
            return len(separated.components) == 1 and separated.components[
                0
            ].radical() == ideal.change_ring(separated.components[0].ring)
    if len(ideal.generators) > 1:
        graph = _graph_prime_decomposition(ideal)
        if graph is not None:
            return True
        if ideal.krull_dimension > 0 and ideal.ring.domain == sp.QQ:
            split = _associated_prime_splitting(ideal, 128)
            if split is not None:
                extended = ideal.change_ring(split.components[0].ring)
                return len(split.components) == 1 and split.components[0] == extended
    if ideal.krull_dimension == 0 and len(ideal.generators) > 1:
        decomposition = absolute_primary_decomposition(ideal)
        return len(decomposition.components) == 1 and decomposition.components[
            0
        ].radical() == ideal.change_ring(decomposition.components[0].ring)
    if len(ideal.generators) == 1 and (
        len(ideal.ring.variables) > 2
        or sp.Poly(ideal.generators[0], *ideal.ring.variables).total_degree() > 2
    ):
        decomposition = absolute_primary_decomposition(ideal)
        return len(decomposition.components) == 1 and decomposition.components[
            0
        ].radical() == ideal.change_ring(decomposition.components[0].ring)
    if len(ideal.generators) != 1 or len(ideal.ring.variables) not in (1, 2):
        from .function_field_constants import algebraic_constant_field, certify_geometric_prime

        try:
            return certify_geometric_prime(ideal)
        except NotImplementedError:
            if ideal.ring.domain != sp.QQ or ideal.krull_dimension <= 0:
                raise
            constants = algebraic_constant_field(ideal)
            return len(constants.geometric_primes) == 1
    polynomial = sp.Poly(ideal.generators[0], *ideal.ring.variables)
    if polynomial.total_degree() == 1:
        return True
    if len(ideal.ring.variables) == 1:
        return False
    if polynomial.total_degree() == 2:
        return _quadric_rank(polynomial, ideal.ring.variables) == 3
    raise NotImplementedError("geometric integrality above degree two is not implemented")


def geometric_components(ideal, *, extension=None):
    """Return certified geometric components of supported finite ideals and hypersurfaces.

    Quadratic splitting fields over QQ are inferred from a discriminant.
    An explicit extension may also be supplied. Every returned factor is
    checked to have degree one, preventing an incomplete decomposition.
    """
    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    if ideal.is_unit:
        return GeometricComponents((), ideal.ring.domain)
    if ideal.is_zero:
        return GeometricComponents((ideal,), ideal.ring.domain)
    if len(ideal.generators) > 1 and ideal.ring.domain == sp.QQ:
        separated = _separated_hypersurface_primary_decomposition(ideal)
        if separated is not None:
            return GeometricComponents(
                tuple(component.radical() for component in separated.components),
                separated.field,
            )
    if len(ideal.generators) > 1:
        graph = _graph_prime_decomposition(ideal)
        if graph is not None:
            return GeometricComponents((ideal,), ideal.ring.domain)
        if ideal.krull_dimension > 0 and ideal.ring.domain == sp.QQ:
            split = _associated_prime_splitting(ideal, 128)
            if split is not None:
                primes = tuple(component.radical() for component in split.components)
                minimal = tuple(
                    prime
                    for prime in primes
                    if not any(other != prime and other.is_subset(prime) for other in primes)
                )
                return GeometricComponents(
                    tuple(
                        prime for index, prime in enumerate(minimal) if prime not in minimal[:index]
                    ),
                    split.field,
                )
    if ideal.krull_dimension == 0 and len(ideal.generators) > 1:
        decomposition = absolute_primary_decomposition(ideal)
        return GeometricComponents(
            tuple(component.radical() for component in decomposition.components),
            decomposition.field,
        )
    if len(ideal.generators) == 1 and (
        len(ideal.ring.variables) > 2
        or sp.Poly(ideal.generators[0], *ideal.ring.variables).total_degree() > 2
    ):
        decomposition = absolute_primary_decomposition(ideal)
        return GeometricComponents(
            tuple(component.radical() for component in decomposition.components),
            decomposition.field,
        )
    if len(ideal.generators) != 1 or len(ideal.ring.variables) not in (1, 2):
        from .function_field_constants import algebraic_constant_field

        constants = algebraic_constant_field(ideal)
        return GeometricComponents(constants.geometric_primes, constants.splitting_field)
    polynomial = sp.Poly(ideal.generators[0], *ideal.ring.variables)
    if len(ideal.ring.variables) == 1 and polynomial.total_degree() > 2:
        decomposition = absolute_primary_decomposition(ideal)
        return GeometricComponents(
            tuple(component.radical() for component in decomposition.components),
            decomposition.field,
        )
    if polynomial.total_degree() > 2:
        raise NotImplementedError(
            "multivariate absolute factor discovery above degree two is not implemented"
        )
    if is_geometrically_integral(ideal):
        return GeometricComponents((ideal.radical(),), ideal.ring.domain)
    if extension is None and ideal.ring.domain == sp.QQ:
        if polynomial.total_degree() == 2 and len(ideal.ring.variables) == 1:
            discriminant = sp.discriminant(polynomial.as_expr(), ideal.ring.variables[0])
            candidate = sp.sqrt(discriminant)
            extension = candidate if candidate.is_Rational is not True else ()
        elif polynomial.total_degree() == 2:
            homogenizer = sp.Dummy("homogenizer")
            homogeneous = polynomial.homogenize(homogenizer).as_expr()
            hessian = sp.hessian(homogeneous, (*ideal.ring.variables, homogenizer)) / 2
            for first in range(3):
                for second in range(first + 1, 3):
                    minor = sp.det(hessian.extract((first, second), (first, second)))
                    if minor != 0:
                        candidate = sp.sqrt(-minor)
                        extension = candidate if candidate.is_Rational is not True else ()
                        break
                if extension is not None:
                    break
        if extension is None:
            extension = ()
    if extension is None:
        raise NotImplementedError(
            "an explicit splitting field is required over this coefficient domain"
        )
    from .ideals import PolynomialRing

    extensions = (extension,) if not isinstance(extension, tuple) else extension
    field = ideal.ring.domain if not extensions else ideal.ring.domain.algebraic_field(*extensions)
    factors = sp.factor_list(
        ideal.generators[0],
        *ideal.ring.variables,
        extension=list(extensions) if extensions else None,
    )[1]
    if any(sp.Poly(f, *ideal.ring.variables).total_degree() != 1 for f, _ in factors):
        raise NotImplementedError("the supplied coefficient extension does not split the quadric")
    ring = PolynomialRing(ideal.ring.variables, field, ideal.ring.order)
    components = tuple(ring.ideal(f) for f, _ in factors)
    return GeometricComponents(components, field)


@dataclass(frozen=True)
class AbsolutePrimaryDecomposition:
    """Geometric primary components after an exact finite algebraic field extension."""

    components: tuple[Ideal, ...]
    field: object

    def verify(self, original, **limits):
        """Independently verify primary status, absolute supports and reconstruction.

        Every associated prime of each component is recomputed over the
        declared splitting field and must equal that component's radical.
        Each radical must have its own geometric-primality certificate.
        Finally, the intersection must exactly recover the extended input,
        and no component may be redundant. Embedded primes are retained.
        No computed certificate, candidate Groebner backend, or component
        classification flag is trusted here.
        """
        from .components import DecompositionLimitError, associated_primes
        from .function_field_constants import certify_geometric_prime
        from .ideals import PolynomialRing

        if not isinstance(original, Ideal):
            return False
        if not self.components:
            return original.is_unit
        ambient = PolynomialRing(original.ring.variables, self.field, original.ring.order)
        expected = original.change_ring(ambient)
        if any(component.ring != ambient or component.is_unit for component in self.components):
            return False
        try:
            primes = []
            for component in self.components:
                prime = component.radical()
                if associated_primes(component, **limits) != (prime,):
                    return False
                if not certify_geometric_prime(prime):
                    return False
                if any(prime == previous for previous in primes):
                    return False
                primes.append(prime)
            result = ambient.ideal(1)
            for component in self.components:
                result = result.intersection(component)
            if result != expected:
                return False
            for index, component in enumerate(self.components):
                others = ambient.ideal(1)
                for offset, other in enumerate(self.components):
                    if index != offset:
                        others = others.intersection(other)
                if others.is_subset(component):
                    return False
            return True
        except (
            DecompositionLimitError,
            NotImplementedError,
            ValueError,
            TypeError,
            ArithmeticError,
        ):
            return False


def _binary_form_primary_decomposition(ideal, max_root_degree):
    """Split homogeneous binary forms into geometric linear primary factors.

    A homogeneous polynomial in two variables becomes a univariate polynomial
    after setting the second coordinate to one. Its geometric factors are
    x-alpha*y and the possibly repeated factor y, with their exact powers.
    """
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or len(ideal.generators) != 1 or len(ideal.ring.variables) != 2:
        return None
    x, y = ideal.ring.variables
    original = sp.Poly(ideal.generators[0], x, y, domain=sp.QQ)
    if not original.is_homogeneous:
        return None
    total_degree = original.total_degree()
    if total_degree > max_root_degree:
        return None
    affine = sp.Poly(original.as_expr().subs(y, 1), x, domain=sp.QQ)
    roots = sp.roots(affine.as_expr(), x)
    if sum(roots.values()) != affine.degree():
        return None
    field = (
        sp.QQ.algebraic_field(*roots)
        if any(sp.sympify(root).is_Rational is not True for root in roots)
        else sp.QQ
    )
    ring = PolynomialRing((x, y), field, ideal.ring.order)
    components = tuple(
        [ring.ideal((x - root * y) ** power) for root, power in roots.items()]
        + (
            [ring.ideal(y ** (total_degree - affine.degree()))]
            if total_degree > affine.degree()
            else []
        )
    )
    reconstructed = ring.ideal(1)
    for component in components:
        reconstructed = reconstructed.intersection(component)
    if reconstructed != ideal.change_ring(ring):
        raise ArithmeticError("homogeneous binary factors failed reconstruction")
    return AbsolutePrimaryDecomposition(components, field)


def _monomial_primary_decomposition(ideal):
    """Split monomial ideals by their coordinate-primary components.

    Monomial primary decomposition is combinatorial and unchanged by extension
    of the coefficient field. Its associated coordinate primes are geometrically
    prime, so this is an absolute decomposition in any dimension.
    """
    from .components import primary_decomposition

    if not ideal.ring.variables or not all(
        len(sp.Poly(g, *ideal.ring.variables).terms()) == 1 for g in ideal.generators
    ):
        return None
    components = primary_decomposition(ideal)
    intersection = ideal.ring.ideal(1)
    for component in components:
        intersection = intersection.intersection(component)
    if intersection != ideal:
        raise ArithmeticError("monomial components failed reconstruction")
    return AbsolutePrimaryDecomposition(components, ideal.ring.domain)


def _projectively_smooth(polynomial, variables):
    """Certify absolute irreducibility using smooth projective completion.

    In projective dimension at least two, distinct hypersurface components
    necessarily intersect. Such an intersection is singular, so smoothness
    of the homogenized projective hypersurface certifies irreducibility over
    the algebraic closure. Checking each affine projective chart excludes the
    irrelevant origin from the Jacobian criterion.
    """
    from .ideals import PolynomialRing

    if len(variables) < 2:
        return False
    original = sp.Poly(polynomial, *variables)
    if original.is_homogeneous and len(variables) >= 3:
        coordinates = variables
        homogeneous = original.as_expr()
    else:
        homogenizer = sp.Dummy("homogenizer")
        homogeneous = original.homogenize(homogenizer).as_expr()
        coordinates = (*variables, homogenizer)
    ambient = PolynomialRing(coordinates, sp.QQ)
    singular = (*[sp.diff(homogeneous, symbol) for symbol in coordinates], homogeneous)
    return all((ambient.ideal(*singular, coordinate - 1)).is_unit for coordinate in coordinates)


def _linear_support_primary_decomposition(ideal):
    """Certify field-relative primary factors with geometrically prime supports.

    A primary ideal stays primary after extending a characteristic-zero
    coefficient field if its associated prime is geometrically prime. An
    affine-linear prime is geometrically prime over every extension, so the
    usual primary decomposition is already absolute in this case.
    """
    from .components import associated_primes, primary_decomposition

    if ideal.ring.domain != sp.QQ:
        return None
    primes = associated_primes(ideal)
    if not primes or any(
        any(sp.Poly(g, *ideal.ring.variables).total_degree() > 1 for g in prime.groebner_basis)
        for prime in primes
    ):
        return None
    components = primary_decomposition(ideal)
    intersection = ideal.ring.ideal(1)
    for component in components:
        intersection = intersection.intersection(component)
    if intersection != ideal:
        raise ArithmeticError("linear-support components failed reconstruction")
    return AbsolutePrimaryDecomposition(components, ideal.ring.domain)


def _algebraic_coordinate_primary_decomposition(ideal, max_root_degree):
    """Split primary factors after discovering algebraic coordinate constants.

    If an associated prime contains a polynomial in one ambient coordinate,
    its roots are algebraic constants on that component. Passing to a common
    number field may split the associated prime. The extended decomposition
    is returned only if *every* associated prime is affine-linear, a simple
    absolute-primality certificate valid over every coefficient extension.
    This includes embedded and nonreduced positive-dimensional components.
    """
    from .components import (
        DecompositionLimitError,
        associated_primes,
        primary_decomposition,
    )
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension <= 0:
        return None
    variables = ideal.ring.variables
    try:
        associated = associated_primes(ideal)
    except (NotImplementedError, DecompositionLimitError):
        return None
    algebraic_constants = []
    for prime in associated:
        for variable in variables:
            eliminated = prime.eliminate(*(other for other in variables if other != variable))
            if eliminated.is_zero:
                continue
            for equation in eliminated.groebner_basis:
                polynomial = sp.Poly(equation, variable, domain=sp.QQ)
                if polynomial.degree() > max_root_degree:
                    return None
                roots = sp.roots(polynomial.as_expr(), variable)
                if sum(roots.values()) != polynomial.degree():
                    return None
                algebraic_constants.extend(root for root in roots if root.is_Rational is not True)
    if not algebraic_constants:
        return None
    if len(algebraic_constants) > 12:
        return None
    field = sp.QQ.algebraic_field(*algebraic_constants)
    ring = PolynomialRing(variables, field, ideal.ring.order)
    extended = ideal.change_ring(ring)
    try:
        primes = associated_primes(extended)
        if not all(
            all(
                sp.Poly(generator, *variables).total_degree() <= 1
                for generator in prime.groebner_basis
            )
            for prime in primes
        ):
            return None
        components = primary_decomposition(extended)
    except (NotImplementedError, DecompositionLimitError):
        return None
    reconstructed = ring.ideal(1)
    for component in components:
        reconstructed = reconstructed.intersection(component)
    if reconstructed != extended:
        raise ArithmeticError("extended primary components failed exact reconstruction")
    return AbsolutePrimaryDecomposition(tuple(components), field)


def _separated_hypersurface_primary_decomposition(ideal, *, max_components=128):
    """Decompose a complete intersection in disjoint variable blocks.

    Absolute factors in disjoint coordinate blocks are geometrically prime
    after taking their sum. Powers of these factors form primary ideals, and
    the distributive intersection reconstructs the complete intersection.
    This is a certificate for a nonprincipal positive-dimensional class; it
    does not claim to handle coupled generators or embedded primes.
    """
    from itertools import product

    from .absolute_factorization import absolute_factor_hypersurface
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or len(ideal.generators) < 2:
        return None
    supports = [set(sp.sympify(generator).free_symbols) for generator in ideal.generators]
    if any(not support for support in supports):
        return None
    if any(left & right for index, left in enumerate(supports) for right in supports[index + 1 :]):
        return None

    blocks = []
    extensions = []
    total = 1
    for generator, support in zip(ideal.generators, supports, strict=True):
        coordinates = tuple(var for var in ideal.ring.variables if var in support)
        if len(coordinates) == 1:
            coordinate = coordinates[0]
            poly = sp.Poly(generator, coordinate, domain=sp.QQ)
            roots = sp.roots(poly.as_expr(), coordinate)
            if sum(roots.values()) != poly.degree():
                return None
            factors = tuple((coordinate - root, power) for root, power in roots.items())
            extensions.extend(root for root in roots if root.is_Rational is not True)
        else:
            factorization = absolute_factor_hypersurface(generator, coordinates)
            factors = factorization.factors
            extensions.extend(getattr(factorization.field, "orig_ext", ()))
        total *= len(factors)
        if total > max_components:
            from .components import DecompositionLimitError

            raise DecompositionLimitError("absolute component count exceeds max_components")
        blocks.append(factors)

    field = sp.QQ.algebraic_field(*extensions) if extensions else sp.QQ
    ring = PolynomialRing(ideal.ring.variables, field, ideal.ring.order)
    components = tuple(
        ring.ideal(*(factor**power for factor, power in choices)) for choices in product(*blocks)
    )
    reconstructed = ring.ideal(1)
    for component in components:
        reconstructed = reconstructed.intersection(component)
    if reconstructed != ideal.change_ring(ring):
        raise ArithmeticError("separated absolute primary decomposition failed reconstruction")
    return AbsolutePrimaryDecomposition(components, field)


def _graph_prime_decomposition(ideal, *, max_candidates=32):
    """Certify geometrically prime polynomial-graph ideals by exact ideal equality.

    Graph relations eliminating each dependent variable identify the
    coordinate ring with a polynomial ring in independent coordinates.
    This criterion is preserved under arbitrary coefficient extensions.
    """
    from itertools import combinations

    from .ideals import PolynomialRing

    variables = ideal.ring.variables
    dimension = ideal.krull_dimension
    if dimension <= 0 or dimension >= len(variables):
        return None
    for index, free_indices in enumerate(combinations(range(len(variables)), dimension)):
        if index >= max_candidates:
            return None
        independent = tuple(variables[pos] for pos in free_indices)
        dependent = tuple(var for var in variables if var not in independent)
        reordered = PolynomialRing((*dependent, *independent), ideal.ring.domain, "lex")
        candidate = ideal.change_ring(reordered)
        if len(candidate.groebner_basis) != len(dependent):
            continue
        representatives = tuple(candidate.normal_form(variable) for variable in dependent)
        if any(value.free_symbols - set(independent) for value in representatives):
            continue
        graph_equations = tuple(
            variable - expression
            for variable, expression in zip(dependent, representatives, strict=True)
        )
        if candidate == reordered.ideal(*graph_equations):
            return AbsolutePrimaryDecomposition((ideal,), ideal.ring.domain)
    return None


def _graph_lift_primary_decomposition(ideal, max_root_degree, max_components):
    """Transport primary components through a certified polynomial graph.

    If the equations have the form (d_i-f_i(u), J(u)), quotienting by
    the graph equations identifies the ambient coordinate ring with
    K[u]. Primary ideals and their associated primes are carried through
    this isomorphism, including embedded components and nilpotents.
    A direct ideal-equality check certifies the proposed presentation.
    """
    from itertools import combinations

    from .ideals import PolynomialRing

    variables = ideal.ring.variables
    if ideal.ring.domain != sp.QQ or len(variables) < 3:
        return None
    attempts = 0
    for size in range(1, len(variables)):
        for chosen in combinations(variables, size):
            attempts += 1
            if attempts > 48:
                return None
            dependent = tuple(variable for variable in variables if variable not in chosen)
            lex = PolynomialRing((*dependent, *chosen), sp.QQ, "lex")
            in_lex = ideal.change_ring(lex)
            representatives = tuple(in_lex.normal_form(variable) for variable in dependent)
            if any(value.free_symbols - set(chosen) for value in representatives):
                continue
            graph_equations = tuple(
                variable - value for variable, value in zip(dependent, representatives, strict=True)
            )
            if not all(in_lex.contains(equation) for equation in graph_equations):
                continue
            substitutions = dict(zip(dependent, representatives, strict=True))
            residual_equations = tuple(
                sp.cancel(equation.subs(substitutions, simultaneous=True))
                for equation in ideal.generators
            )
            base_ring = PolynomialRing(chosen, sp.QQ, ideal.ring.order)
            base = base_ring.ideal(*residual_equations)
            if base.is_zero or base.is_unit:
                continue
            if in_lex != lex.ideal(*graph_equations, *base.generators):
                continue
            try:
                result = absolute_primary_decomposition(
                    base, max_root_degree=max_root_degree, max_components=max_components
                )
            except (NotImplementedError, ValueError):
                continue
            extended = PolynomialRing(variables, result.field, ideal.ring.order)
            components = tuple(
                extended.ideal(*graph_equations, *component.generators)
                for component in result.components
            )
            recovered = extended.ideal(1)
            for component in components:
                recovered = recovered.intersection(component)
            if recovered != ideal.change_ring(extended):
                raise ArithmeticError("graph-lifted primary components failed reconstruction")
            return AbsolutePrimaryDecomposition(components, result.field)
    return None


def _associated_prime_splitting(ideal, max_components):
    """Split associated primes after discovering an algebraic constant extension.

    Polynomial factors in associated-prime Gröbner bases provide possible
    splitting constants. After extending the *whole* ideal, redo primary
    decomposition, retaining embedded components. Accept only if every
    extended associated prime has a geometric-primality certificate, each
    primary component has the expected associated prime, and all components
    reconstruct the extended ideal.
    """
    from .absolute_factorization import FactorizationLimitError, absolute_factor_hypersurface
    from .components import DecompositionLimitError, associated_primes, primary_decomposition
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension <= 0:
        return None
    variables = ideal.ring.variables
    try:
        associated = associated_primes(ideal)
    except (DecompositionLimitError, NotImplementedError):
        return None
    roots = []
    for prime in associated:
        for expression in prime.groebner_basis:
            polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
            if polynomial.total_degree() <= 1 or polynomial.total_degree() > 6:
                continue
            try:
                factors = absolute_factor_hypersurface(
                    expression, variables, max_unknowns=80, max_roots=12
                )
            except (FactorizationLimitError, NotImplementedError, ValueError):
                continue
            roots.extend(getattr(factors.field, "orig_ext", ()))
            if len(roots) > 8:
                return None
    if not roots:
        return None
    field = sp.QQ.algebraic_field(*roots)
    ring = PolynomialRing(variables, field, ideal.ring.order)
    extended = ideal.change_ring(ring)
    try:
        primary = primary_decomposition(extended)
        primes = associated_primes(extended)
    except (DecompositionLimitError, NotImplementedError):
        return None
    if len(primary) > max_components:
        raise DecompositionLimitError("absolute component count exceeds max_components")
    if not all(_certified_geometric_prime(prime) for prime in primes):
        return None
    for component in primary:
        # This also excludes components that are merely primary to a
        # different field-relative prime after the coefficient extension.
        prime = component.radical()
        if prime not in primes or associated_primes(component) != (prime,):
            raise ArithmeticError("extended primary component failed associated-prime check")
    recovered = ring.ideal(1)
    for component in primary:
        recovered = recovered.intersection(component)
    if recovered != extended:
        raise ArithmeticError("extended primary components failed reconstruction")
    return AbsolutePrimaryDecomposition(tuple(primary), field)


def _certified_geometric_prime(ideal):
    """Check geometric primality through exact sufficient certificates."""
    from .function_field_constants import certify_geometric_prime

    try:
        return certify_geometric_prime(ideal)
    except NotImplementedError:
        return False


def _constant_field_primary_decomposition(ideal, max_components):
    """Split nonprincipal ideals using proven full generic constant fields.

    The construction preserves embedded components by extending the original
    ideal and recomputing all primary factors over the common splitting field.
    Every resulting associated prime must have its own geometric-primality
    certificate; equality of the intersection alone is insufficient.
    """
    from .components import DecompositionLimitError, associated_primes, primary_decomposition
    from .function_field_constants import algebraic_constant_field
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension <= 0:
        return None
    try:
        primes = associated_primes(ideal)
    except (DecompositionLimitError, NotImplementedError):
        return None
    roots = []
    for prime in primes:
        if prime.krull_dimension <= 0:
            continue
        try:
            constants = algebraic_constant_field(prime, max_candidates=12)
        except (NotImplementedError, ValueError):
            return None
        roots.extend(getattr(constants.splitting_field, "orig_ext", ()))
    if not roots:
        # Do not add a second, potentially expensive primary-decomposition
        # path when no coefficient extension is required.
        return None
    field = sp.QQ.algebraic_field(*roots)
    ambient = PolynomialRing(ideal.ring.variables, field, ideal.ring.order)
    extended = ideal.change_ring(ambient)
    try:
        components = primary_decomposition(extended)
        associated = associated_primes(extended)
    except (DecompositionLimitError, NotImplementedError):
        return None
    if len(components) > max_components:
        raise DecompositionLimitError("absolute component count exceeds max_components")
    if not all(_certified_geometric_prime(prime) for prime in associated):
        return None
    reconstructed = ambient.ideal(1)
    for component in components:
        prime = component.radical()
        if prime not in associated or associated_primes(component) != (prime,):
            raise ArithmeticError("a primary component has an unexpected associated prime")
        reconstructed = reconstructed.intersection(component)
    if reconstructed != extended:
        raise ArithmeticError("constant-field primary factors failed reconstruction")
    return AbsolutePrimaryDecomposition(tuple(components), field)


def _smooth_hypersurface_primary_decomposition(ideal):
    """Split products of powers of projectively smooth hypersurfaces."""
    if ideal.ring.domain != sp.QQ or len(ideal.generators) != 1 or len(ideal.ring.variables) < 2:
        return None
    variables = ideal.ring.variables
    factors = sp.Poly(ideal.generators[0], *variables, domain=sp.QQ).factor_list()[1]
    if not all(
        sp.Poly(factor, *variables).total_degree() == 1
        or _projectively_smooth(factor.as_expr(), variables)
        for factor, _ in factors
    ):
        return None
    components = tuple(
        ideal.ring.ideal(factor.as_expr() ** exponent) for factor, exponent in factors
    )
    intersection = ideal.ring.ideal(1)
    for component in components:
        intersection = intersection.intersection(component)
    if intersection != ideal:
        raise ArithmeticError("hypersurface components failed reconstruction")
    return AbsolutePrimaryDecomposition(components, ideal.ring.domain)


def absolute_primary_decomposition(ideal, *, max_root_degree=4, max_components=128):
    """Compute exact absolute primary factors when geometric supports can be certified.

    A univariate polynomial over QQ splits into powers of linear maximal ideals
    over its constructed splitting field; all multiplicities are preserved.
    Principal ideals use exact absolute hypersurface factorization. Certain
    nonprincipal ideals are handled through disjoint-variable complete
    intersections, monomial decomposition, and certified geometric supports.
    The max_components bound limits the number of separated intersections.
    General nonprincipal absolute primary decomposition remains unsupported.
    """
    from .ideals import PolynomialRing

    if not isinstance(ideal, Ideal):
        raise TypeError("an Ideal is required")
    max_components = integer(max_components, "max_components", 1)
    if ideal.is_unit:
        return AbsolutePrimaryDecomposition((), ideal.ring.domain)
    if ideal.is_zero:
        return AbsolutePrimaryDecomposition((ideal,), ideal.ring.domain)
    if ideal.ring.domain == sp.QQ and len(ideal.generators) == 1 and len(ideal.ring.variables) >= 2:
        from .absolute_factorization import (
            FactorizationLimitError,
            absolute_factor_hypersurface,
        )

        try:
            factorization = absolute_factor_hypersurface(
                ideal.generators[0],
                ideal.ring.variables,
                max_unknowns=100,
                max_roots=max(12, max_root_degree),
            )
        except FactorizationLimitError:
            # The smooth-projective criterion can still establish geometric
            # primality when the residue matrix is beyond its degree budget.
            smooth = _smooth_hypersurface_primary_decomposition(ideal)
            if smooth is not None:
                return smooth
            raise
        ambient = PolynomialRing(ideal.ring.variables, factorization.field, ideal.ring.order)
        components = tuple(ambient.ideal(factor**power) for factor, power in factorization.factors)
        intersection = ambient.ideal(1)
        for component in components:
            intersection = intersection.intersection(component)
        if intersection != ideal.change_ring(ambient):
            raise ArithmeticError("absolute factors failed primary ideal reconstruction")
        return AbsolutePrimaryDecomposition(components, factorization.field)
    binary = _binary_form_primary_decomposition(ideal, max_root_degree)
    if binary is not None:
        return binary
    monomial = _monomial_primary_decomposition(ideal)
    if monomial is not None:
        return monomial
    separated = _separated_hypersurface_primary_decomposition(ideal, max_components=max_components)
    if separated is not None:
        return separated
    if len(ideal.generators) > 1:
        graph = _graph_prime_decomposition(ideal)
        if graph is not None:
            return graph
        lifted = _graph_lift_primary_decomposition(ideal, max_root_degree, max_components)
        if lifted is not None:
            return lifted
    smooth = _smooth_hypersurface_primary_decomposition(ideal)
    if smooth is not None:
        return smooth
    # Higher-dimensional schemes supported on affine-linear geometric primes
    # have no hidden absolute splitting after a coefficient extension.
    if len(ideal.ring.variables) >= 2 and ideal.krull_dimension > 0:
        algebraic = _algebraic_coordinate_primary_decomposition(ideal, max_root_degree)
        if algebraic is not None:
            return algebraic
        linear = _linear_support_primary_decomposition(ideal)
        if linear is not None:
            return linear
        associated_split = _associated_prime_splitting(ideal, max_components)
        if associated_split is not None:
            return associated_split
        constant_split = _constant_field_primary_decomposition(ideal, max_components)
        if constant_split is not None:
            return constant_split
    if len(ideal.generators) != 1:
        try:
            return _shape_primary_decomposition(ideal, max_root_degree)
        except NotImplementedError:
            return _artinian_primary_decomposition(ideal, max_root_degree)
    variables = ideal.ring.variables
    if len(variables) == 1 and ideal.ring.domain == sp.QQ:
        polynomial = sp.Poly(ideal.generators[0], variables[0], domain=sp.QQ)
        if polynomial.degree() > max_root_degree:
            raise NotImplementedError("polynomial degree exceeds max_root_degree")
        roots = sp.roots(polynomial.as_expr(), variables[0])
        if sum(roots.values()) != polynomial.degree():
            raise NotImplementedError(
                "complete exact radical representation of roots is unavailable"
            )
        field = sp.QQ.algebraic_field(*roots) if roots else sp.QQ
        ring = PolynomialRing(variables, field, ideal.ring.order)
        components = tuple(
            ring.ideal((variables[0] - root) ** power) for root, power in roots.items()
        )
    elif len(variables) <= 2 and sp.Poly(ideal.generators[0], *variables).total_degree() <= 2:
        geometric = geometric_components(ideal)
        field = geometric.field
        ring = PolynomialRing(variables, field, ideal.ring.order)
        factors = sp.factor_list(
            ideal.generators[0],
            *variables,
            extension=list(getattr(field, "orig_ext", ()))
            if getattr(field, "is_AlgebraicField", False)
            else None,
        )[1]
        components = tuple(ring.ideal(factor**power) for factor, power in factors)
    else:
        raise NotImplementedError(
            "absolute primary decomposition is available only for low-degree hypersurfaces or bounded univariate polynomials"
        )
    expected = ideal.change_ring(ring)
    intersection = ring.ideal(1)
    for component in components:
        intersection = intersection.intersection(component)
    if intersection != expected:
        raise ArithmeticError("absolute components did not reconstruct the input ideal")
    return AbsolutePrimaryDecomposition(components, field)


def _shape_primary_decomposition(ideal, max_root_degree):
    """Split finite monogenic quotient algebras presented in lexicographic shape."""
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension != 0 or not ideal.ring.variables:
        raise NotImplementedError(
            "absolute primary decomposition requires a finite shape-position ideal over QQ"
        )
    variables = ideal.ring.variables
    last = variables[-1]
    first = variables[:-1]
    lex = ideal.change_ring(PolynomialRing(variables, sp.QQ, "lex"))
    representatives = tuple(lex.normal_form(v) for v in first)
    if any(rep.free_symbols - {last} for rep in representatives):
        raise NotImplementedError("finite ideal does not have the declared primitive-element shape")
    minimal = lex.eliminate(*first)
    if len(minimal.groebner_basis) != 1:
        raise NotImplementedError("the primitive element has no principal minimal polynomial")
    polynomial = sp.Poly(minimal.groebner_basis[0], last, domain=sp.QQ)
    if polynomial.degree() > max_root_degree:
        raise NotImplementedError("primitive element degree exceeds max_root_degree")
    roots = sp.roots(polynomial.as_expr(), last)
    if sum(roots.values()) != polynomial.degree():
        raise NotImplementedError("complete exact root representation is unavailable")
    field = sp.QQ.algebraic_field(*roots)
    ring = PolynomialRing(variables, field, ideal.ring.order)
    components = tuple(
        ring.ideal(
            *(v - representative for v, representative in zip(first, representatives, strict=True)),
            (last - root) ** power,
        )
        for root, power in roots.items()
    )
    recovered = ring.ideal(1)
    for component in components:
        recovered = recovered.intersection(component)
    if recovered != ideal.change_ring(ring):
        raise ArithmeticError("shape-position primary components do not reconstruct the ideal")
    return AbsolutePrimaryDecomposition(components, field)


def _primitive_primary_decomposition(ideal, max_root_degree):
    """Find a separating primitive element for a finite monogenic quotient algebra.

    A fresh coordinate t=linear(x) converts a zero-dimensional quotient
    into shape position if the selected linear form generates the algebra.
    Each conjugate root then determines an exact primary component, obtained
    by eliminating t. In particular, this supports radical finite ideals
    without a shape-position coordinate in the original presentation.
    """
    from .ideals import PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension != 0:
        raise NotImplementedError("primitive-element decomposition requires a finite ideal over QQ")
    variables = ideal.ring.variables
    parameter = sp.Dummy("primitive")
    extended = PolynomialRing((*variables, parameter), sp.QQ, "lex")
    candidates = (
        *(variable for variable in reversed(variables)),
        *(
            sum((index + 1) ** power * variable for index, variable in enumerate(variables))
            for power in range(1, 7)
        ),
    )
    for candidate in candidates:
        graph = extended.ideal(*ideal.generators, parameter - candidate)
        representatives = tuple(graph.normal_form(variable) for variable in variables)
        if any(rep.free_symbols - {parameter} for rep in representatives):
            continue
        elimination = graph.eliminate(*variables)
        if len(elimination.groebner_basis) != 1:
            continue
        polynomial = sp.Poly(elimination.groebner_basis[0], parameter, domain=sp.QQ)
        if polynomial.degree() > max_root_degree:
            raise NotImplementedError("primitive element degree exceeds max_root_degree")
        roots = sp.roots(polynomial.as_expr(), parameter)
        if sum(roots.values()) != polynomial.degree():
            raise NotImplementedError("complete algebraic root representation is unavailable")
        field = sp.QQ.algebraic_field(*roots)
        ambient = PolynomialRing(variables, field, ideal.ring.order)
        lifted = PolynomialRing((*variables, parameter), field, "lex")
        components = tuple(
            lifted.ideal(
                *(
                    variable - representative
                    for variable, representative in zip(variables, representatives, strict=True)
                ),
                (parameter - root) ** power,
            )
            .eliminate(parameter)
            .change_ring(ambient)
            for root, power in roots.items()
        )
        recovered = ambient.ideal(1)
        for component in components:
            recovered = recovered.intersection(component)
        if recovered != ideal.change_ring(ambient):
            raise ArithmeticError("primitive-element primary factors do not reconstruct the ideal")
        return AbsolutePrimaryDecomposition(components, field)
    raise NotImplementedError(
        "the finite algebra is not monogenic under the searched primitive elements"
    )


def _artinian_primary_decomposition(ideal, max_root_degree):
    """Split a finite algebra by exact eigenvalues of a separating linear form.

    The multiplication operator of a linear form has characteristic polynomial
    with roots equal to its values at geometric points, repeated according to
    local algebra lengths. Distinct eigenvalues separate primary factors even
    if the local factors are not monogenic. A factor is certified local by
    checking that its radical has geometric residue degree one over the
    constructed splitting field.
    """
    from .ideals import CoordinateRing, PolynomialRing

    if ideal.ring.domain != sp.QQ or ideal.krull_dimension != 0:
        raise NotImplementedError(
            "finite primary splitting requires a zero-dimensional ideal over QQ"
        )
    algebra = CoordinateRing(ideal).finite_algebra(max_dimension=max_root_degree)
    variables = ideal.ring.variables
    parameter = sp.Symbol("_alggeom_spectrum")
    if parameter in variables:
        raise ValueError("reserved spectrum symbol conflicts with a polynomial coordinate")
    candidates = (
        *reversed(variables),
        *(
            sum((index + 1) ** power * variable for index, variable in enumerate(variables))
            for power in range(1, 9)
        ),
    )
    for linear in candidates:
        matrix = algebra.multiplication_matrix(linear)
        characteristic = sp.Poly(sp.Matrix(matrix).charpoly(parameter).as_expr(), parameter)
        roots = sp.roots(characteristic.as_expr(), parameter)
        if sum(roots.values()) != algebra.dimension:
            continue
        field = sp.QQ if all(root.is_Rational for root in roots) else sp.QQ.algebraic_field(*roots)
        ring = PolynomialRing(variables, field, ideal.ring.order)
        components = tuple(
            ring.ideal(*ideal.generators, (linear - root) ** power) for root, power in roots.items()
        )
        if any(
            CoordinateRing(component.radical())
            .finite_algebra(max_dimension=max_root_degree)
            .dimension
            != 1
            for component in components
        ):
            continue
        recovered = ring.ideal(1)
        for component in components:
            recovered = recovered.intersection(component)
        if recovered != ideal.change_ring(ring):
            raise ArithmeticError("Artinian primary factors do not reconstruct the ideal")
        return AbsolutePrimaryDecomposition(components, field)
    raise NotImplementedError("no separating linear form was found within the search bound")


__all__ = [
    "GeometricComponents",
    "AbsolutePrimaryDecomposition",
    "is_geometrically_integral",
    "geometric_components",
    "absolute_primary_decomposition",
]
