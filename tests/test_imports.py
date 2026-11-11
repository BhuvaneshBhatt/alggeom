"""The core package remains usable without optional real geometry imports."""

import subprocess
import sys


def test_core_imports(tmp_path):
    program = f"""
import importlib.abc
import sys
sys.path[:] = {sys.path!r}
class CoreOnly(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] == 'semialg' or fullname in {{
            'algroots.functional_decomposition', 'algroots._algebraic_discovery',
            'algroots._multivariate_decomposition', 'algroots._decomposition_chains'
        }}:
            raise ImportError('unexpected dependency: '+fullname)
sys.meta_path.insert(0, CoreOnly())
import sympy as sp
from alggeom.decomposition import decompose, decompose_relation
x,z=sp.symbols('x z')
maps=decompose(x**8,x)
value=x
for component in reversed(maps):
    value=component(value)
assert sp.cancel(value-x**8) == 0
assert len(decompose_relation(x**8-z,(x,z))) == 3
assert 'semialg' not in sys.modules
"""
    process = subprocess.run(
        [sys.executable, "-I", "-B", "-c", program],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert process.returncode == 0, process.stdout + process.stderr
