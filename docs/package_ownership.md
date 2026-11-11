# Package ownership

Use the package that owns the mathematical question. The three packages share exact arithmetic and keep their dependency direction acyclic.

| Question | Canonical package | Main interfaces |
| --- | --- | --- |
| Gröbner-basis construction, finite quotient arithmetic, exact coefficient fields, zero-dimensional solving | [algroots](https://github.com/BhuvaneshBhatt/algroots) | `quotient`, `domains`, `finite_fields`, `finite_decomposition`, root solvers |
| Reusable coefficient-field and plane-curve arithmetic | [algroots](https://github.com/BhuvaneshBhatt/algroots) | `subfields`, `function_fields`, `curves` |
| General coordinate rings, schemes, normalization and function decomposition | [alggeom](https://github.com/BhuvaneshBhatt/alggeom) | `ideals`, `rings`, `schemes`, `normalization`, `decomposition` |
| Real loci, inequalities, parameter conditions and exact real graph verification | [semialg](https://github.com/BhuvaneshBhatt/semialg) | `decision`, `qe`, `conditional`, `function_graph`, `function_verification` |

algroots has no dependency on alggeom or semialg. semialg depends on algroots; it interprets exact algebra through real signs, inequalities and witnesses. alggeom depends on algroots and uses semialg through its optional real extra. Cross-package applications belong in alggeom's real profile. The backends do not import the application package.

A positive-dimensional coordinate ring has a Krull dimension and generally infinite dimension over its coefficient field. `algroots.QuotientAlgebra.dimension` is a finite vector-space dimension. A real branch condition belongs to a semialgebraic graph, not to an algebraic ideal or its Zariski closure. Finite fields do not have a real ordering.

Function-decomposition search, coefficient reconstruction, chains and their branch-aware adapters are supplied by `alggeom.decomposition`. The backend field arithmetic is reusable independently of a decomposition search. Generic real graph existence, uniqueness and composition checks remain available through `semialg.function_verification`.

The real decomposition adapters call semialg's generic graph-comparison engine and return its verification records. Rational-coordinate membership and coefficient rebasing use the public `algroots.rational_coordinates` interface. The same finite quotient backend supplies the arithmetic used by semialg's GTZ certificates.

Alggeom has no native F4 kernel. Its ideal and quotient factories call algroots and retain the requested backend policy. The independent `groebner_certificates` module checks basis and membership evidence for geometric use; parameter-cover algorithms organize these calculations over geometric strata. Kernel arithmetic, scheduling and matrix elimination are tested in algroots. Adapter routing and geometric proof obligations are tested in alggeom.
