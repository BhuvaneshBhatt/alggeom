"""Compare geometric primary components with Macaulay2's independent engine.

Set M2_EXECUTABLE to the Macaulay2 executable if it is not on PATH. Set
M2_PREFIX when running a manually extracted Debian distribution. Macaulay2
itself must be installed separately; ordinary test runs skip these cases.
"""

import os
import shutil
import signal
import subprocess
import tempfile
from contextlib import suppress
from types import SimpleNamespace

import pytest
import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing

x, y, z = sp.symbols("x y z")


def _m2_polynomial(expression):
    """Serialize fixture polynomials in Macaulay2 arithmetic notation."""
    return sp.sstr(sp.expand(expression)).replace("**", "^").replace("I", "a")


def _m2_ideal(expressions):
    return f"ideal({','.join(_m2_polynomial(expr) for expr in expressions)})"


def _m2_command():
    executable = os.environ.get("M2_EXECUTABLE") or shutil.which("M2")
    if not executable:
        pytest.skip("Macaulay2 executable not available")
    command = [executable]
    if prefix := os.environ.get("M2_PREFIX"):
        command.extend(("--prefix", prefix))
    command.extend(("--no-debug", "--no-readline", "--no-prompts", "--no-tty", "-q"))
    return command


def _run_m2(script):
    """Run an oracle in its own process group and clean up its descendants."""
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as output:
        process = subprocess.Popen(
            _m2_command(),
            stdin=subprocess.PIPE,
            stdout=output,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        try:
            process.communicate(input=script, timeout=60)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.communicate()
            raise
        finally:
            # Macaulay2 may leave subprocesses holding inherited descriptors.
            with suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGTERM)
        output.seek(0)
        return SimpleNamespace(returncode=process.returncode, stdout=output.read())


def test_macaulay2_detects_geometric_splitting():
    """A QQ-prime need not stay prime after algebraic extension."""
    common = "I=ideal(x^2+y^2,z)\n"
    checks = (
        ("R=QQ[x,y,z]\n" + common, "true"),
        ("K=QQ[a]/(a^2+1)\nR=K[x,y,z]\n" + common, "false"),
    )
    for setup, expected in checks:
        run = _run_m2(setup + 'print("ALGGEOM_GEOMETRIC_PRIME=" | toString (isPrime I))\nexit 0\n')
        assert run.returncode == 0, run.stdout
        assert f"ALGGEOM_GEOMETRIC_PRIME={expected}" in run.stdout
        assert ": error:" not in run.stdout


@pytest.mark.parametrize(
    ("name", "generators", "geometric_extension"),
    [
        ("embedded_double_line", (x**2, x * y), False),
        ("nonprincipal_components", (x * y, x * z), False),
        ("mixed_components", (x**2, x * y, y * z), False),
        ("complex_embedded", (x**2 + 1, y**2, y * z), True),
        ("geometric_nonprincipal", (x**2 + y**2, z), True),
        pytest.param(
            "geometric_nonreduced",
            ((x**2 + y**2) ** 2, z),
            True,
            marks=pytest.mark.skipif(
                os.environ.get("ALGGEOM_SLOW_ORACLE") != "1",
                reason="expensive algebraic-extension oracle; opt in with ALGGEOM_SLOW_ORACLE=1",
            ),
        ),
    ],
)
def test_macaulay2_primary_components(name, generators, geometric_extension):
    """Check primariness, associated primes and reconstruction independently."""
    ring = PolynomialRing((x, y, z))
    decomposition = absolute_primary_decomposition(ring.ideal(*generators))
    components = tuple(component.generators for component in decomposition.components)
    # These fixtures have affine-linear radical supports. This makes primality
    # over the tested number field sufficient for geometric primality; it would
    # not be sufficient for a general nonlinear prime ideal.
    assert all(
        all(
            sp.Poly(g, x, y, z, extension=True).total_degree() <= 1
            for g in component.radical().groebner_basis
        )
        for component in decomposition.components
    )
    coeff_field = "K=QQ[a]/(a^2+1)\nR=K[x,y,z]" if geometric_extension else "R=QQ[x,y,z]"
    script = (
        "\n".join(
            (
                coeff_field,
                f"I={_m2_ideal(generators)}",
                f"C={{{','.join(_m2_ideal(c) for c in components)}}}",
                "P=primaryDecomposition I",
                'print("ALGGEOM_ORACLE_PRIMARY=" | toString (all(C, J -> isPrimary J)))',
                'print("ALGGEOM_ORACLE_GEOMETRIC=" | toString (all(C, J -> isPrime radical J)))',
                'print("ALGGEOM_ORACLE_RECONSTRUCT=" | toString (intersect C == I))',
                'print("ALGGEOM_ORACLE_ASSOCIATED=" | toString '
                "(#C == #P and all(associatedPrimes I, Q -> "
                "any(C, J -> radical(J) == Q))))",
                'print("ALGGEOM_ORACLE_COMPONENTS=" | toString '
                "(#C == #P and all(C, J -> any(P, Q -> J == Q))))",
                "exit 0",
            )
        )
        + "\n"
    )
    run = _run_m2(script)
    required = (
        "ALGGEOM_ORACLE_PRIMARY=true",
        "ALGGEOM_ORACLE_GEOMETRIC=true",
        "ALGGEOM_ORACLE_RECONSTRUCT=true",
        "ALGGEOM_ORACLE_ASSOCIATED=true",
        "ALGGEOM_ORACLE_COMPONENTS=true",
    )
    assert run.returncode == 0, f"Macaulay2 failed for {name}:\n{run.stdout}"
    assert all(text in run.stdout for text in required), (
        f"Macaulay2 disagrees with alggeom for {name}:\n{run.stdout}"
    )
    assert ": error:" not in run.stdout, f"Macaulay2 reported an error for {name}:\n{run.stdout}"
