"""Keep generated references and executable benchmark reports reproducible."""

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_tool(path, arguments, directory):
    """Run a tool outside the project while retaining the current installation under test."""
    environment = {
        **os.environ,
        "PYTHONPATH": os.pathsep.join(str(Path(p).resolve()) for p in sys.path),
    }
    return subprocess.run(
        [sys.executable, "-B", str(ROOT / path), *arguments],
        cwd=directory,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )


@pytest.mark.real
def test_api_reference(tmp_path):
    result = run_tool("tools/generate_api.py", ["--check"], tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_benchmark_report(tmp_path):
    destination = tmp_path / "measurements.json"
    result = run_tool(
        "benchmarks/run.py",
        [
            "--case",
            "factor_discovery",
            "--suite",
            "smoke",
            "--repeat",
            "2",
            "--warmup",
            "0",
            "--output",
            str(destination),
        ],
        tmp_path,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(destination.read_text(encoding="utf-8"))
    assert report["repeat"] == 2 and report["warmup"] == 0 and report["omitted_cases"] == []
    assert report["environment"]["python"] == sys.version.split()[0]
    assert set(report["environment"]["source_sha256"]) == {"alggeom", "algroots"}
    (case,) = report["benchmarks"]
    assert case["status"] == "ok" and case["name"] == "factor_discovery"
    assert case["dimensions"] == {
        "output_degree": 3,
        "factor_degree": 2,
        "factor_algebra_dimension": 3,
    }
    assert len(case["samples_seconds"]) == 2
    assert all(sample >= 0 for sample in case["samples_seconds"])
    assert case["median_seconds"] == sum(case["samples_seconds"]) / 2


def test_backend_inputs():
    inputs = tomllib.loads((ROOT / "ci/backends.toml").read_text(encoding="utf-8"))
    assert inputs == {
        "algroots": {"requirement": "algroots>=0.4.0"},
        "semialg": {"requirement": "semialg>=1.5.0"},
    }
