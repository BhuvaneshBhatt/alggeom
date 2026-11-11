# Getting started

Use alggeom when equations and maps carry geometric structure: a multiple fiber, a singular local ring, a finite projection or a normalization. The example below studies a cusp using its defining ideal and local invariants.

## Install the package

Use Python 3.11 or later:

```sh
python -m pip install alggeom
```

Install `"alggeom[real]"` when a calculation also needs inequalities, real domains or selected branches. This extra requires semialg 1.5.0 or newer; the core package requires algroots 0.4.0 or newer.

Gröbner construction uses the native algroots backend by default. Compatible wheels include its compiled kernels. The [backend guide](backends.md) explains native availability and explicit SymPy selection.

## Study a singular curve

<!-- examples: executable -->
```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme

x, y = sp.symbols("x y")
S = PolynomialRing((x, y))
C = AffineScheme(S.ideal(y**2 - x**3))
assert C.tangent_cone((0, 0)) == S.ideal(y**2)
assert C.local_dimension((0, 0)) == 1
assert C.embedding_dimension((0, 0)) == 2
assert C.coordinate_ring.vector_space_dimension is sp.oo
```

The cusp is a curve, so its local dimension is one. Its two-dimensional tangent space detects the singularity at the origin, and its coordinate ring has infinite vector-space dimension over the coefficient field. Finite quotient dimension and geometric dimension answer different questions.

## Choose the next calculation

Continue with [worked applications](worked_applications.md) for fibers, normalization and real branches. Read [conventions](conventions.md) before comparing dimensions, closures or component counts, and check the [coefficient-field matrix](coefficient_fields_support.md) before changing fields. The [documentation index](index.md) groups the mathematical guides; the [API reference](api.md) gives exact signatures and result types.
