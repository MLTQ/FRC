"""Existing execution budgets around the released-workspace full workflow."""
import argparse
from . import guarded_tensor_refinement as guard
from .released_metric_rates import run as calculate
from .implicit_plasma_records import provenance


def run(destination, control, nr, nz, original=None):
    saved = guard.calculate, guard.provenance
    try:
        guard.calculate = calculate
        guard.provenance = lambda seeds: provenance([*seeds, 'examples.guarded_released_refinement',
            'examples.guarded_tensor_refinement', 'examples.released_metric_rates'])
        return guard.run(destination, control, nr, nz, original)
    finally:
        guard.calculate, guard.provenance = saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination'); parser.add_argument('control')
    parser.add_argument('--nr', type=int, required=True); parser.add_argument('--nz', type=int, required=True)
    parser.add_argument('--original'); args = parser.parse_args()
    run(args.destination, args.control, args.nr, args.nz, args.original)
