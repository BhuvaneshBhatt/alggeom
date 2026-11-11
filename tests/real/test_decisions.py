"""Mathematical Boolean entry points requiring real geometry."""

import sympy as sp

from alggeom.decomposition import decomposition_condition, real_equivalent
from alggeom.decomposition.chains import rational_chain

x, y, a = sp.symbols("x y a", real=True)


def test_boolean_results():
    condition = decomposition_condition(x**4 + a * x, x, degrees=(2, 2), parameters=(a,))
    assert sp.simplify(sp.Equivalent(condition, sp.Eq(a, 0))) is sp.true
    assert real_equivalent(rational_chain(x**8, x), output=y) is True
    assert (
        real_equivalent(
            rational_chain(x**8, x), output=y, target_selector=y < 0, check_functions=False
        )
        is False
    )
