"""Reuse the unchanged physical experiment with the coordinate-addition solver."""
import argparse
import json
from pathlib import Path
from plasma_bench.lifted_modal_gauss import LiftedModalGauss
from . import spitzer_bessel_decay as experiment


def run(destination, *, nr=8, steps=80, amplitude=1e-4):
    original_solver = experiment.GuardedModalGauss
    original_provenance = experiment.provenance
    try:
        experiment.GuardedModalGauss = LiftedModalGauss
        experiment.provenance = lambda seeds: original_provenance([*seeds, 'examples.lifted_spitzer_decay'])
        record = experiment.run(destination, nr=nr, steps=steps, amplitude=amplitude)
    finally:
        experiment.GuardedModalGauss = original_solver
        experiment.provenance = original_provenance
    record['integrator_variant'] = 'LiftedModalGauss: coordinate addition before magnetic lift; original equations and bounds'
    (Path(destination)/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--nr', type=int, default=8); parser.add_argument('--steps', type=int, default=80)
    parser.add_argument('--amplitude', type=float, default=1e-4)
    args = parser.parse_args(); run(args.destination, nr=args.nr, steps=args.steps, amplitude=args.amplitude)
