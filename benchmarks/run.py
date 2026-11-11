"""Measure decomposition workloads in separate processes and write reproducible JSON."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

from cases import CASES, REAL_CASES, prepare
from sympy.external.gmpy import GROUND_TYPES


def module_digest(name):
    """Fingerprint a package or benchmark directory, including its relative file paths."""
    root = Path(importlib.import_module(name).__file__).parent
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py"), key=lambda path: path.relative_to(root).as_posix()):
        digest.update(path.relative_to(root).as_posix().encode() + b"\0")
        digest.update(path.read_bytes().replace(b"\r\n", b"\n") + b"\0")
    return digest.hexdigest()


def environment(include_real):
    """Record interpreter, platform, arithmetic backend and installed source fingerprints."""
    names = ["alggeom", "algroots", "sympy", "python-flint"]
    if include_real:
        names.append("semialg")
    versions = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "ground_types": GROUND_TYPES,
        "benchmark_sha256": module_digest("cases"),
        "packages": versions,
        "source_sha256": {
            name: module_digest(name)
            for name in ("alggeom", "algroots", *(("semialg",) if include_real else ()))
        },
    }


def measure(name, suite, repeat, warmup):
    """Time only operations, checking every result after the clock has stopped."""
    started = time.perf_counter()
    workload = prepare(name, suite)
    setup = time.perf_counter() - started
    samples = []
    for index in range(warmup + repeat):
        started = time.perf_counter()
        result = workload.operation()
        elapsed = time.perf_counter() - started
        if not workload.validate(result):
            raise RuntimeError(f"{name}: benchmark result failed its mathematical check")
        if index >= warmup:
            samples.append(elapsed)
    return {
        "name": name,
        "status": "ok",
        "dimensions": workload.dimensions,
        "setup_seconds": setup,
        "samples_seconds": samples,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
    }


def main():
    """Run selected cases with per-case timeouts; failed workloads never become timing samples."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("smoke", "standard"), default="standard")
    parser.add_argument("--case", choices=CASES, action="append", dest="cases")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument(
        "--timeout",
        type=float,
        default=120,
        help="per-case process limit, including setup and checks",
    )
    parser.add_argument("--real", action="store_true", help="include the optional real workload")
    parser.add_argument("--output", type=Path, default=Path("benchmark.json"))
    parser.add_argument("--worker", choices=CASES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repeat < 1 or args.warmup < 0 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("repeat and timeout must be positive; warmup must be nonnegative")
    if args.worker:
        print(json.dumps(measure(args.worker, args.suite, args.repeat, args.warmup)))
        return
    if args.real and importlib.util.find_spec("semialg") is None:
        parser.error("--real requires semialg; install the matching real backend")
    selected = args.cases or list(CASES)
    omitted = [name for name in selected if name in REAL_CASES and not args.real]
    if args.cases and omitted:
        parser.error("the real_branches case requires --real")
    cases = [name for name in selected if name not in omitted]
    records = []
    for name in cases:
        command = [
            sys.executable,
            "-S",
            "-B",
            str(Path(__file__).resolve()),
            "--worker",
            name,
            "--suite",
            args.suite,
            "--repeat",
            str(args.repeat),
            "--warmup",
            str(args.warmup),
        ]
        try:
            process = subprocess.run(
                command,
                env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)},
                capture_output=True,
                text=True,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            record = {"name": name, "status": "timeout", "timeout_seconds": args.timeout}
        else:
            record = (
                json.loads(process.stdout)
                if process.returncode == 0
                else {"name": name, "status": "failed", "error": process.stderr[-4000:]}
            )
        records.append(record)
        print(f"{name}: {record['status']}", flush=True)
    report = {
        "recorded_at": datetime.now(UTC).isoformat(),
        "environment": environment(args.real),
        "suite": args.suite,
        "repeat": args.repeat,
        "warmup": args.warmup,
        "cache_policy": "one fresh process per case; repeated operations retain normal SymPy and backend caches",
        "omitted_cases": omitted,
        "benchmarks": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    if any(record["status"] != "ok" for record in records):
        parser.exit(1, f"Benchmark failure; details saved to {args.output}\n")


if __name__ == "__main__":
    main()
