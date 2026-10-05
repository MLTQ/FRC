"""Immutable measured steps for the Newton-product continuation solver."""
from . import matched_step_benchmark as benchmark
from plasma_bench.recycled_modal_gauss import RecycledModalGauss
from .implicit_plasma_records import provenance


def run(destination, spec, *, epsilon=1e-4, intervals=60, steps=3):
    original_class = benchmark.AcceleratedModalGauss
    original_provenance = benchmark.provenance
    try:
        benchmark.AcceleratedModalGauss = RecycledModalGauss
        benchmark.provenance = lambda seeds: provenance([*seeds, 'examples.recycled_step_benchmark'])
        return benchmark.run(destination, spec, epsilon=epsilon, intervals=intervals, steps=steps)
    finally:
        benchmark.AcceleratedModalGauss = original_class
        benchmark.provenance = original_provenance
