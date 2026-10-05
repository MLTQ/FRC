"""Immutable measured steps for the checked continuation solver."""
from . import matched_step_benchmark as benchmark
from plasma_bench.continued_modal_gauss import ContinuedModalGauss
from .implicit_plasma_records import provenance


def run(destination, spec, *, epsilon=1e-4, intervals=60, steps=3):
    # The underlying protocol's seed is extended explicitly for this variant.
    original_class = benchmark.AcceleratedModalGauss
    original_provenance = benchmark.provenance
    try:
        benchmark.AcceleratedModalGauss = ContinuedModalGauss
        benchmark.provenance = lambda seeds: provenance([*seeds, 'examples.continued_step_benchmark'])
        return benchmark.run(destination, spec, epsilon=epsilon, intervals=intervals, steps=steps)
    finally:
        benchmark.AcceleratedModalGauss = original_class
        benchmark.provenance = original_provenance
