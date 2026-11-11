"""Recover a hidden scalar substitution from an expanded rational map."""

import sympy as sp

from alggeom.decomposition import decompose
from alggeom.decomposition.rational import search_rational

x, y, t = sp.symbols("x y t")
inner = (x * x + y * y) / (x * y + 1)
function = sp.cancel((inner * inner + 1) / (inner * inner + inner + 2))
search = search_rational(function, (x, y), parameter=t)
assert search.complete and search.verify()
assert len(search.decompositions) == 1
pair = search.decompositions[0]
assert pair.degrees == (2, 2)
assert sp.cancel(pair.outer.subs(t, pair.inner) - function) == 0
maps = decompose((x * x + y * y) ** 4, (x, y))
assert len(maps) == 3
print("Outer-to-inner maps:", maps)
