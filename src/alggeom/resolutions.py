"""Minimal standard-graded polynomial-module resolutions and Betti invariants."""

from collections import Counter
from dataclasses import dataclass

import sympy as sp

from ._algebra import integer
from .graded import _parameter
from .ideals import Ideal, PolynomialRing
from .modules import PresentedModule, Submodule, syzygies


def _degrees(ring, matrix, row_degrees):
    """Infer homogeneous column shifts, including the shifts of the target free module."""
    degrees = []
    for j in range(matrix.cols):
        terms = {
            sum(monomial) + row_degrees[i]
            for i in range(matrix.rows)
            for monomial, _ in ring._poly(matrix[i, j]).items()
        }
        if len(terms) > 1:
            raise ValueError("the presentation must be homogeneous with its generator degrees")
        degrees.append(next(iter(terms)) if terms else None)
    return tuple(degrees)


def _minimal_columns(ring, matrix, row_degrees):
    """Choose a homogeneous minimal generating set by graded membership.

    Processing increasing degrees makes membership in earlier columns the
    graded Nakayama test. Equal-degree dependencies have scalar coefficients.
    """
    degrees = _degrees(ring, matrix, row_degrees)
    selected = sp.ImmutableMatrix(matrix.rows, 0, [])
    shifts = []
    indices = sorted((j for j, d in enumerate(degrees) if d is not None), key=degrees.__getitem__)
    for j in indices:
        if not Submodule(ring, selected).contains(matrix[:, j]):
            selected = selected.row_join(matrix[:, j])
            shifts.append(degrees[j])
    return selected, tuple(shifts)


def _minimal_presentation(module, generator_degrees):
    """Cancel scalar pivots before minimizing the homogeneous relation columns.

    A pivot relation expresses its row generator in the other generators.
    The Schur complement substitutes that expression into all other relations;
    the retained original generators give the augmentation to the input module.
    """
    ring = module.ring
    matrix = module.presentation
    shifts = list(generator_degrees)
    augmentation = sp.ImmutableMatrix.eye(matrix.rows)
    _degrees(ring, matrix, shifts)
    while True:
        pivot = next(
            (
                (i, j)
                for i in range(matrix.rows)
                for j in range(matrix.cols)
                if matrix[i, j] != 0 and not any(map(sum, ring._poly(matrix[i, j]).monoms()))
            ),
            None,
        )
        if pivot is None:
            break
        i, j = pivot
        rows = [k for k in range(matrix.rows) if k != i]
        columns = [k for k in range(matrix.cols) if k != j]
        scalar = matrix[i, j]
        matrix = sp.ImmutableMatrix(
            len(rows),
            len(columns),
            [
                ring.convert(matrix[row, column] - matrix[row, j] * matrix[i, column] / scalar)
                for row in rows
                for column in columns
            ],
        )
        augmentation = augmentation[:, rows]
        shifts.pop(i)
    matrix, relation_degrees = _minimal_columns(ring, matrix, tuple(shifts))
    return matrix, tuple(shifts), relation_degrees, augmentation


@dataclass(frozen=True)
class GradedResolution:
    """A minimal graded complex F_i -> F_(i-1) augmenting to a presented module.

    ``shifts[i]`` lists j for the summands S(-j) in F_i. ``differentials[i-1]``
    maps F_i to F_(i-1). The augmentation maps F_0 to the original presentation
    generators. A truncated complex has ``complete=False`` and does not supply
    projective dimension, regularity or a module Hilbert series.
    """

    module: PresentedModule
    generator_degrees: tuple[int, ...]
    differentials: tuple[sp.ImmutableMatrix, ...]
    shifts: tuple[tuple[int, ...], ...]
    augmentation: sp.ImmutableMatrix
    complete: bool

    @property
    def betti_numbers(self):
        """The dictionary beta_(i,j), indexed by homological and internal degrees."""
        return dict(Counter((i, j) for i, degrees in enumerate(self.shifts) for j in degrees))

    @property
    def betti_table(self):
        """Rows j-i mapped to tuples of beta_(i,j), with homological degree in columns.

        A truncated complex supplies the known Betti columns through its last
        term; higher homological degrees have not been computed.
        """
        betti = self.betti_numbers
        return {
            row: tuple(betti.get((i, i + row), 0) for i in range(len(self.shifts)))
            for row in sorted({j - i for i, j in betti})
        }

    @property
    def projective_dimension(self):
        """The module's projective dimension; the zero module has value -1."""
        self._require_complete()
        return max((i for i, shifts in enumerate(self.shifts) if shifts), default=-1)

    @property
    def regularity(self):
        """Castelnuovo-Mumford regularity; the zero module has value -oo."""
        self._require_complete()
        return max((j - i for i, shifts in enumerate(self.shifts) for j in shifts), default=-sp.oo)

    def hilbert_series(self, *, parameter=None):
        """Return the alternating graded free-module Hilbert series of the resolved module."""
        self._require_complete()
        parameter = _parameter(self.module.ring.ideal(), parameter, "t")
        numerator = sum(
            (-1) ** i * parameter**j for i, shifts in enumerate(self.shifts) for j in shifts
        )
        return sp.cancel(numerator / (1 - parameter) ** len(self.module.ring.variables))

    def _require_complete(self):
        if not self.complete:
            raise ValueError("this invariant requires a complete resolution")

    def verify(self):
        """Check grading, minimality, augmentation and exactness by independent kernel tests.

        For a truncated complex, exactness is checked at every term except
        its leftmost term. Completeness additionally requires a zero terminal
        kernel. Verification recomputes syzygies and can cost as much as construction.
        """
        ring = self.module.ring
        if not isinstance(ring, PolynomialRing) or len(self.shifts) != len(self.differentials) + 1:
            return False
        if not self.shifts or self.augmentation.shape != (
            self.module.generator_count,
            len(self.shifts[0]),
        ):
            return False
        if len(self.generator_degrees) != self.module.generator_count:
            return False
        try:
            _degrees(ring, self.module.presentation, self.generator_degrees)
            augmentation_degrees = _degrees(ring, self.augmentation, self.generator_degrees)
        except ValueError:
            return False
        if augmentation_degrees != self.shifts[0]:
            return False
        for i, matrix in enumerate(self.differentials, 1):
            if matrix.shape != (len(self.shifts[i - 1]), len(self.shifts[i])):
                return False
            try:
                degrees = _degrees(ring, matrix, self.shifts[i - 1])
            except ValueError:
                return False
            if degrees != self.shifts[i]:
                return False
            if any(
                value != 0 and not any(map(sum, ring._poly(value).monoms())) for value in matrix
            ):
                return False
        presentation = self.module.presentation
        augmented = Submodule(ring, presentation.row_join(self.augmentation))
        if not all(
            augmented.contains(sp.eye(presentation.rows)[:, j]) for j in range(presentation.rows)
        ):
            return False
        kernel = syzygies(self.augmentation.row_join(-presentation), ring)[: len(self.shifts[0]), :]
        previous = self.augmentation
        for i, matrix in enumerate(self.differentials):
            image = Submodule(ring, matrix)
            if not all(image.contains(kernel[:, j]) for j in range(kernel.cols)):
                return False
            product = previous * matrix
            if i == 0:
                if not all(
                    self.module.relations.contains(product[:, j]) for j in range(product.cols)
                ):
                    return False
            elif any(ring.convert(v) != 0 for v in product):
                return False
            kernel = syzygies(matrix, ring)
            previous = matrix
        return not self.complete or all(v == 0 for v in kernel)


def free_resolution(module, *, generator_degrees=None, max_length=None):
    """Compute a minimal standard-graded free resolution over S=K[x_1,...,x_n].

    An Ideal means the quotient S/I. A PresentedModule must be over S; use
    ``restrict_scalars()`` explicitly for a module over S/I. Integer generator
    shifts default to zero. Relations must be homogeneous for these shifts.
    Successive exact syzygies are minimized by graded membership. The default
    length bound is n, the Hilbert syzygy bound. A smaller explicit bound
    returns a truncated complex marked incomplete when another term is needed.
    """
    if isinstance(module, Ideal):
        module = PresentedModule(
            module.ring, sp.ImmutableMatrix(1, len(module.generators), module.generators)
        )
    if not isinstance(module, PresentedModule):
        raise TypeError("an Ideal or PresentedModule is required")
    if not isinstance(module.ring, PolynomialRing):
        raise ValueError("resolve over a polynomial ring; call restrict_scalars() for S/I modules")
    generator_degrees = (
        (0,) * module.generator_count if generator_degrees is None else tuple(generator_degrees)
    )
    if len(generator_degrees) != module.generator_count:
        raise ValueError("one generator degree is required per presentation row")
    degrees = tuple(integer(d, "generator degree", minimum=None) for d in generator_degrees)
    bound = len(module.ring.variables) if max_length is None else integer(max_length, "max_length")
    matrix, degrees, relation_degrees, augmentation = _minimal_presentation(module, degrees)
    matrices = []
    shifts = [degrees]
    while matrix.cols and len(matrices) < bound:
        matrices.append(matrix)
        shifts.append(relation_degrees)
        kernel = syzygies(matrix, module.ring)
        matrix, relation_degrees = _minimal_columns(module.ring, kernel, relation_degrees)
    complete = not matrix.cols
    if max_length is None and not complete:
        raise RuntimeError("the graded syzygies exceeded the polynomial-ring length bound")
    return GradedResolution(
        module, tuple(generator_degrees), tuple(matrices), tuple(shifts), augmentation, complete
    )


__all__ = ["GradedResolution", "free_resolution"]
