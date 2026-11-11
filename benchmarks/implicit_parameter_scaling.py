"""Repeated exact algebra benchmarks for implicit cubic and quartic parameters.

Run using ``python benchmarks/implicit_parameter_scaling.py`` with alggeom
installed. Timing values depend on the arithmetic backend and host; neither
correctness tests nor CI compare raw durations across different machines.
"""

import json
import platform
import statistics
from time import perf_counter

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.parameter_cover import comprehensive_prime_cover
from alggeom.parameter_fields import PrimeParameterField, groebner_over_prime_parameters


def measure(call, repeats=5):
    call()  # Cache initialization and interpreter warmup are excluded.
    observations = []
    for _ in range(repeats):
        start = perf_counter()
        result = call()
        observations.append(perf_counter() - start)
        if result is None:
            raise ArithmeticError("benchmark returned no algebraic result")
    return {"median_seconds": statistics.median(observations), "samples": observations}


def main():
    a, b, x, y = sp.symbols("a b x y")
    output = {
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "host": platform.platform(),
        "workloads": {},
    }
    for degree in (2, 3, 4):
        equation = a**degree + b**3 + 1
        parameter_ideal = PolynomialRing((a, b), sp.QQ).ideal(equation)

        def arithmetic(parameter_ideal=parameter_ideal):
            field = PrimeParameterField(parameter_ideal)
            value = field.element(a**2 + 1, b + 1)
            if value * value.inverse() != field.element(1):
                raise ArithmeticError("quotient field identity failed")
            return value

        def generic_basis(parameter_ideal=parameter_ideal):
            certificate = groebner_over_prime_parameters(
                (a * x + b * y, y - 1), (x, y), parameter_ideal
            )
            if not certificate.verify():
                raise ArithmeticError("generic quotient basis certificate failed")
            return certificate

        def recursive_cover(equation=equation):
            cover = comprehensive_prime_cover((equation * x, a * y), (x, y), (a, b), max_nodes=25)
            if not cover.verify():
                raise ArithmeticError("parameter-cover certificate failed")
            return cover

        for label, call in (
            ("quotient_arithmetic", arithmetic),
            ("quotient_groebner", generic_basis),
            ("recursive_cover", recursive_cover),
        ):
            output["workloads"][f"degree_{degree}_{label}"] = measure(call)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
