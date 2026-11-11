# Read search results correctly

Choose an interface according to the answer you need. A composition identity, an exhaustive search result, and equality of real functions answer different questions.

| Question | Interface | What to inspect |
| --- | --- | --- |
| Give one fully resolved rational chain | `decompose` | Lambda tuple, or an incomplete-result exception |
| Give every rational two-map class | `decompositions` | Pair tuple, or an incomplete-result exception |
| Keep partial rational evidence | `search_rational` | `decompositions`, `examined`, `candidate_count`, `complete` |
| Examine implicit coefficient-field decompositions | `search_algebraic` | Factor coverage, reconstruction attempts, stated search scope |
| Keep a recursive proof tree | `rational_chain`, `algebraic_chain` | `leaves`, `stop_reason`, `complete`, `verify()` |
| Preserve a selected real function | `verify_decomposition`, `verify_chain` | `equivalent`, totality, uniqueness, domain and selectors |

## Finding a split versus finishing a search

A successful candidate proves existence as soon as its exact identity verifies. Other candidates may remain unexamined. For `x**6`, a budget of one finds a split while leaving other two-map classes unresolved.

```python
import sympy as sp
from alggeom.decomposition.rational import search_rational

x = sp.Symbol("x")
partial = search_rational(x**6, x, max_candidates=1)
assert partial.status == "decomposable"
assert not partial.complete
assert partial.examined == 1 and partial.candidate_count == 3
assert partial.verify()
full = search_rational(x**6, x)
assert full.complete and len(full.decompositions) == 2
```

An empty decomposition tuple proves a negative result only when coverage is complete in the stated scope. In particular, an algebraic search's scope is narrower than arbitrary algebraic correspondence decomposition. Consult the [capability limits](capabilities.md) before treating any status as a global mathematical conclusion.

## Work with an incomplete chain

The mathematical entry points require the complete result they promise. When a limit prevents that result, `IncompleteDecompositionError.result` retains the exact partial certificate. It can still be useful for substitutions or real graph checks.

```python
import sympy as sp
from alggeom.decomposition import IncompleteDecompositionError, decompose
from alggeom.decomposition.chains import rational_chain

x = sp.Symbol("x")
try:
    decompose(x**8, x, max_depth=1)
except IncompleteDecompositionError as error:
    partial = error.result
    assert partial.verify() and not partial.complete
    assert partial.compose() == x**8
    assert {leaf.stop_reason for leaf in partial.leaves} == {"max_depth"}
else:
    raise AssertionError("the depth limit must leave unresolved leaves")
resolved = rational_chain(x**8, x, max_depth=4)
assert resolved.complete and len(resolved.components) == 3
```

Increasing a limit restarts the search with a larger budget; the API does not currently resume a saved cursor. A complete selected chain also does not enumerate all possible chains. Finally, field composition can cancel an internal pole: use the real verification interfaces when preserving the affine domain matters.
