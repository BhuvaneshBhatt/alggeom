# Graded resolutions and Betti tables

A free resolution describes all relations among module generators, relations among those relations, and so on. Its graded Betti numbers count the generators required in each homological and polynomial degree. These invariants help compare ideals that have the same dimension and reveal the structure behind Hilbert series.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.resolutions import free_resolution

x, y, z, w, t = sp.symbols("x y z w t")
S = PolynomialRing((x, y, z, w))
twisted_cubic = S.ideal(x * z - y * y, x * w - y * z, y * w - z * z)
resolution = free_resolution(twisted_cubic)
assert resolution.shifts == ((0,), (2, 2, 2), (3, 3))
assert resolution.betti_numbers == {(0, 0): 1, (1, 2): 3, (2, 3): 2}
assert resolution.betti_table == {0: (1, 0, 0), 1: (0, 3, 2)}
assert resolution.projective_dimension == 2
assert resolution.regularity == 1
assert sp.cancel(resolution.hilbert_series(parameter=t) - (1 + 2 * t) / (1 - t) ** 2) == 0
assert resolution.verify()
```

An `Ideal` argument means the module S/I. For a `PresentedModule`, rows are generators and columns are relations. `generator_degrees` supplies the shifts of the original generators, defaulting to zero. A summand S(-j) has its generator in degree j. `shifts[i]` lists those degrees in homological degree i, and `differentials[i-1]` maps that free module to the preceding one. `augmentation` maps its degree-zero free term to the original presentation generators.

`betti_numbers[(i,j)]` is beta_(i,j). `betti_table` uses rows j-i and columns i; absent entries are zero. `regularity` is the maximum j-i. The zero module has projective dimension -1 and regularity -oo. The module's projective dimension counts resolution length; it is different from `graded.projective_dimension(I)`, the geometric dimension of Proj(S/I).

## Base rings and grading

```python
from alggeom.modules import PresentedModule

P = PolynomialRing((x,))
dual_numbers = PresentedModule(P.quotient(x * x), sp.zeros(1, 0))
as_polynomial_module = dual_numbers.restrict_scalars()
finite = as_polynomial_module.free_resolution()
assert finite.shifts == ((0,), (2,))
assert finite.hilbert_series(parameter=t) == 1 + t
```

Resolutions are minimal for the standard polynomial grading over S=K[x_1,...,x_n]. Integer generator shifts, including negative shifts, are supported; relations must be homogeneous for those shifts. Scalar pivots cancel redundant presentation generators, homogeneous membership removes redundant relations, and successive exact syzygies supply the remaining differentials. Hilbert's syzygy theorem bounds the length by n.

Modules over S/I must explicitly call `restrict_scalars()` to obtain an S-presentation. A resolution over S/I may be infinite; this implementation does not substitute a finite S-resolution for that different object. Weighted and multigraded resolutions, general Ext and Tor functors, sheaf cohomology, and resolutions over general quotient rings are not supplied. The [local intersection interface](local_geometry.md) separately computes the Koszul Tor groups required by its supported complete-intersection case.

## Verification and truncation

`verify()` recomputes module kernels, image containments, chain identities, augmentation, grading and minimality. A chain of matrices whose consecutive products vanish need not be exact; verification checks the stronger kernel-image condition. It can cost as much as constructing the resolution.

An explicit `max_length` can return a prefix with `complete=False`. Its known Betti columns are available, and verification excludes exactness at the leftmost uncomputed kernel. Projective dimension, regularity and a module Hilbert series require a complete resolution and raise on a prefix. These are exact algorithms for small polynomial presentations; syzygies can grow substantially before the length bound is reached. See [Eisenbud's commutative-algebra text](https://link.springer.com/book/10.1007/978-1-4612-5350-1) for the algebraic conventions.
