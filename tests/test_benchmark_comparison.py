"""Separate performance-regression decisions from mathematical correctness tests."""

import importlib.util
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[1] / "benchmarks" / "compare.py"
spec = importlib.util.spec_from_file_location("benchmark_compare", MODULE)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def _report(times, *, python="3.13.5", status="ok", dimensions=None):
    dimensions = dimensions or {"parameters": 2}
    return {
        "suite": "smoke",
        "environment": {
            "python": python,
            "implementation": "CPython",
            "machine": "x86_64",
            "ground_types": "flint",
            "packages": {
                "sympy": "1.14.0",
                "python-flint": "0.9.0",
                "algroots": "0.3.0",
            },
        },
        "benchmarks": [
            {
                "name": "groebner_nested",
                "status": status,
                "dimensions": dimensions,
                "samples_seconds": times,
                "median_seconds": sorted(times)[len(times) // 2],
            }
        ],
    }


def test_regression_comparison():
    baseline = _report([1.0, 1.0, 1.0])
    faster = _report([0.7, 0.8, 0.8])
    slower = _report([1.7, 1.8, 1.8])
    assert not compare.compare_reports(baseline, faster)[0]["regression"]
    assert compare.compare_reports(baseline, slower)[0]["regression"]


def test_regression_refuses_incompatible_measurements():
    baseline = _report([1.0, 1.0, 1.0])
    for other in (
        _report([1.0, 1.0, 1.0], python="3.14.0"),
        _report([1.0, 1.0, 1.0], dimensions={"parameters": 3}),
        _report([1.0, 1.0, 1.0], status="timeout"),
        _report([1.0]),
    ):
        with pytest.raises(compare.ComparisonError):
            compare.compare_reports(baseline, other)
