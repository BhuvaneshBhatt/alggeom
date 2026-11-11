# Mathematical conventions

Coordinate rings and maps retain nilpotents. `CoordinateRing.krull_dimension` measures geometric dimension; `vector_space_dimension` is infinite for a nonempty positive-dimensional affine algebra. A nonempty projective scheme has dimension one less than its homogeneous coordinate ring. The empty scheme has dimension -1 here.

Graph ideals are saturated before elimination; image ideals describe Zariski closures, which may contain missing image points. A selected real branch and an input pole guard are extra semialgebraic data.

Noether normalization finds a finite polynomial subring. Integral closure adds missing integral functions and supplies a finite birational map. Plane-curve normalization analysis computes projective invariants and rationality data; an affine integral-closure presentation is a different result.

A module presentation has generators in rows and relations in columns. Kernels, images and cokernels have explicit maps, but their chosen generators need not be unique. A complex is a resolution only when image equals kernel; vanishing consecutive products alone is insufficient.

See [coefficient fields](coefficient_fields_support.md), [geometry](geometry.md), [local geometry](local_geometry.md), [modules](modules.md), [normalization](normalization.md) and [capabilities](capabilities.md) for the exact supported cases.
