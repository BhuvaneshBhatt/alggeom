"""Compare paired benchmark reports on the same execution environment."""

import argparse
import json
import math
from pathlib import Path


class ComparisonError(ValueError):
    """Benchmark reports cannot be compared as paired measurements."""


def compare_reports(baseline, candidate, *, max_ratio=1.5):
    """Return measured ratios; reject incomparable environments or failed workloads.

    Source fingerprints may differ because this tool compares revisions. Machine,
    backend, interpreter and major dependency versions must match instead.
    Each case must provide at least two recorded timings to reduce sensitivity
    to single-run timing noise.
    """
    if not math.isfinite(max_ratio) or max_ratio <= 0:
        raise ValueError("max_ratio must be positive and finite")
    if baseline["suite"] != candidate["suite"]:
        raise ComparisonError("benchmark suites differ")
    for field in ("python", "implementation", "machine", "ground_types"):
        if baseline["environment"].get(field) != candidate["environment"].get(field):
            raise ComparisonError(f"benchmark environments differ: {field}")
    for package in ("sympy", "python-flint", "algroots"):
        old = baseline["environment"]["packages"].get(package)
        new = candidate["environment"]["packages"].get(package)
        if old != new:
            raise ComparisonError(f"benchmark dependencies differ: {package}")
    old_cases = {item["name"]: item for item in baseline["benchmarks"]}
    new_cases = {item["name"]: item for item in candidate["benchmarks"]}
    if old_cases.keys() != new_cases.keys() or not old_cases:
        raise ComparisonError("benchmark case sets differ or are empty")
    changes = []
    for name in sorted(old_cases):
        old, new = old_cases[name], new_cases[name]
        if old["status"] != "ok" or new["status"] != "ok":
            raise ComparisonError(f"failed or incomplete benchmark: {name}")
        if old["dimensions"] != new["dimensions"]:
            raise ComparisonError(f"structural input differs: {name}")
        if len(old["samples_seconds"]) < 2 or len(new["samples_seconds"]) < 2:
            raise ComparisonError(f"insufficient timing samples: {name}")
        earlier, later = old["median_seconds"], new["median_seconds"]
        if earlier <= 0 or later <= 0:
            raise ComparisonError(f"invalid median: {name}")
        ratio = later / earlier
        changes.append(
            {
                "name": name,
                "baseline": earlier,
                "candidate": later,
                "ratio": ratio,
                "regression": ratio > max_ratio,
            }
        )
    return changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--max-ratio", type=float, default=1.5)
    args = parser.parse_args()
    old = json.loads(args.baseline.read_text(encoding="utf-8"))
    new = json.loads(args.candidate.read_text(encoding="utf-8"))
    try:
        changes = compare_reports(old, new, max_ratio=args.max_ratio)
    except ComparisonError as exc:
        parser.exit(2, str(exc) + "\n")
    for change in changes:
        print(
            f"{change['name']}: {change['ratio']:.2f}x "
            f"({'review' if change['regression'] else 'within limit'})"
        )
    if any(change["regression"] for change in changes):
        parser.exit(1, "Performance regression review required\n")


if __name__ == "__main__":
    main()
