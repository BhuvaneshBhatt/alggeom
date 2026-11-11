# Examples

Every catalog entry is an executable application. The installed-example check runs each in a fresh interpreter outside the checkout and executes marked tutorial blocks.

| Script | Problem | Profile |
| --- | --- | --- |
| [affine_normalization.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/affine_normalization.py) | affine normalization | core |
| [algebraic_relations.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/algebraic_relations.py) | algebraic relations | core |
| [coefficient_fields.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/coefficient_fields.py) | coefficient fields | core |
| [curve_geometry.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/curve_geometry.py) | curve geometry | core |
| [curve_normalization.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/curve_normalization.py) | curve normalization | core |
| [curve_reparametrization.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/curve_reparametrization.py) | curve reparametrization | core |
| [cusp_branch.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/cusp_branch.py) | cusp branch | real |
| [finite_points.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/finite_points.py) | finite points | core |
| [graded_resolution.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/graded_resolution.py) | graded resolution | core |
| [linkage_workspace.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/linkage_workspace.py) | linkage workspace | real |
| [local_geometry.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/local_geometry.py) | local geometry | core |
| [module_arithmetic.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/module_arithmetic.py) | module arithmetic | core |
| [module_rank.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/module_rank.py) | module rank | core |
| [number_fields.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/number_fields.py) | number fields | core |
| [parameter_families.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/parameter_families.py) | parameter families | real |
| [primary_components.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/primary_components.py) | primary components | core |
| [quadratic_fibers.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/quadratic_fibers.py) | quadratic fibers | real |
| [quadratic_towers.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/quadratic_towers.py) | quadratic stages | real |
| [rational_maps.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/rational_maps.py) | rational maps | core |
| [real_branches.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/real_branches.py) | real branches | real |
| [response_contours.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/response_contours.py) | response contours | real |
| [scheme_fibers.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/scheme_fibers.py) | scheme fibers | core |
| [solve_equation.py](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/solve_equation.py) | solve equation | core |

Run `python tools/check_examples.py --profile real` after installing the wheel. Core alggeom workflows can use `--profile core`. Every listed script contains assertions for its stated result.
