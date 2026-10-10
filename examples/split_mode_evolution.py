"""Complete original Gauss trajectory using opt-in split-inverse induction."""
from . import matched_mode_evolution as evolution
from .implicit_plasma_records import provenance
from plasma_bench.guarded_modal_gauss import GuardedModalGauss
from plasma_bench.split_inverse_wall_plasma import split_inverse


def run(destination, spec, **options):
    original = evolution.RecycledModalGauss, evolution.preconditioned, evolution.provenance
    try:
        evolution.RecycledModalGauss = GuardedModalGauss
        evolution.preconditioned = split_inverse
        evolution.provenance = lambda seeds: provenance([*seeds, 'examples.split_mode_evolution'])
        return evolution.run(destination, spec, **options)
    finally:
        evolution.RecycledModalGauss, evolution.preconditioned, evolution.provenance = original
