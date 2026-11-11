"""Degree, finiteness and critical loci of exact affine polynomial maps."""

import sympy as sp

from .ideals import CoordinateRing, PolynomialRing
from .schemes import AffineMorphism


def is_dominant(mapping):
    """Test density of an affine morphism's image in its reduced target."""
    if not isinstance(mapping, AffineMorphism):
        raise TypeError("an affine morphism is required")
    return mapping.image_closure.ideal.radical() == mapping.target.ideal.radical()


def _polynomial_graph(mapping):
    if not isinstance(mapping, AffineMorphism):
        raise TypeError("an affine morphism is required")
    if mapping.source.domain != mapping.target.domain:
        raise ValueError("morphism invariants require one common coefficient field")
    ring = PolynomialRing(
        (*mapping.source.variables, *mapping.target.variables), mapping.source.domain, "lex"
    )
    return ring.ideal(
        *mapping.source.ideal.generators,
        *(
            coordinate - expression
            for coordinate, expression in zip(
                mapping.target.variables, mapping.coordinates, strict=True
            )
        ),
    )


def is_finite(mapping):
    """Test module finiteness of an affine morphism by its graph initial ideal.

    With source variables preceding target variables in lexicographic order,
    the graph coordinate ring is finite over the target coordinate ring exactly
    when its initial ideal contains a pure power of every source variable.
    The criterion also applies to nonreduced targets and positive characteristic.
    """
    graph = _polynomial_graph(mapping)
    if mapping.source.is_empty:
        return True
    source_count = len(mapping.source.variables)
    leading = (basis.LM for basis in graph._basis)
    powers = tuple(tuple(monomial) for monomial in leading)
    return all(
        any(
            exponents[index] > 0
            and all(power == 0 for j, power in enumerate(exponents) if j != index)
            for exponents in powers
        )
        for index in range(source_count)
    )


def _generic_fiber(mapping, *, max_dimension=4096):
    """Return the exact zero-dimensional generic fiber over the target function field."""
    if not isinstance(mapping, AffineMorphism):
        raise TypeError("an affine morphism is required")
    if not mapping.target.ideal.is_zero or mapping.target.domain.characteristic() != 0:
        raise NotImplementedError(
            "generic degree currently requires an affine-space target in characteristic zero"
        )
    if not is_dominant(mapping):
        raise ValueError("morphism has no generically finite dominant degree")
    coordinates = mapping.source.variables
    target = mapping.target.variables
    domain = mapping.source.domain
    if target:
        domain = domain.frac_field(*target)
    ring = PolynomialRing(coordinates, domain)
    generic = ring.ideal(
        *mapping.source.ideal.generators,
        *(
            expression - coordinate
            for expression, coordinate in zip(mapping.coordinates, target, strict=True)
        ),
    )
    if generic.krull_dimension != 0:
        raise ValueError("morphism is not generically finite")
    algebra = CoordinateRing(generic).finite_algebra(max_dimension=max_dimension)
    return int(algebra.dimension)


def morphism_degree(mapping, *, max_dimension=4096):
    """Return the generic finite-fiber length of a dominant affine morphism.

    The target must be affine space over a characteristic-zero exact field.
    Generic length agrees with geometric fiber cardinality when the generic
    fiber is reduced, but counts multiplicities for nonreduced fibers.
    The morphism need not be finite: generically finite open maps have degree.
    """
    return _generic_fiber(mapping, max_dimension=max_dimension)


def is_birational(mapping):
    """Test whether an integral dominant morphism induces degree-one fields."""
    if not mapping.source.is_integral or not mapping.target.is_integral:
        raise ValueError("birationality requires integral source and target")
    if not is_dominant(mapping):
        return False
    return morphism_degree(mapping) == 1


def ramification_locus(mapping):
    """Return the non-unramified locus of a finite affine morphism.

    The closed subscheme is defined by Fitt_0(Omega_{B/A}), where B is the
    source coordinate algebra and A is the target coordinate algebra. Its
    presentation has one row per source coordinate and columns given by the
    Jacobian of the source equations and the pullbacks of target coordinates.
    This also applies when the source or target is singular or nonreduced.
    It detects failure of unramifiedness, not the more general failure of
    smoothness or flatness.
    """
    _polynomial_graph(mapping)
    if not is_finite(mapping):
        raise ValueError("ramification requires a finite affine morphism")
    variables = mapping.source.variables
    source_equations = mapping.source.ideal.generators
    relations = (*source_equations, *mapping.coordinates)
    matrix = sp.Matrix(
        len(variables),
        len(relations),
        lambda row, column: sp.diff(relations[column], variables[row]),
    )
    from .modules import PresentedModule, fitting_ideal

    module = PresentedModule(mapping.source.coordinate_ring, matrix)
    return mapping.source.subscheme(*fitting_ideal(module, 0).generators)


def branch_locus(mapping):
    """Return the closed image of the non-unramified locus of a finite affine map."""
    critical = ramification_locus(mapping)
    return critical.morphism(mapping.target, mapping.coordinates).image_closure


__all__ = [
    "is_dominant",
    "is_finite",
    "morphism_degree",
    "is_birational",
    "ramification_locus",
    "branch_locus",
]
