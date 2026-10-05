"""Complete accepted-step protocol with bounded proposal work."""
import argparse
import json
from . import matched_mode_evolution as evolution
from plasma_bench.guarded_modal_gauss import GuardedModalGauss
from .implicit_plasma_records import provenance


def run(destination, spec, **options):
    original_class = evolution.RecycledModalGauss; original_provenance = evolution.provenance
    try:
        evolution.RecycledModalGauss = GuardedModalGauss
        evolution.provenance = lambda seeds: provenance([*seeds, 'examples.guarded_mode_evolution'])
        return evolution.run(destination, spec, **options)
    finally:
        evolution.RecycledModalGauss = original_class; evolution.provenance = original_provenance


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('spec')
    parser.add_argument('--epsilon',type=float,default=1e-4);parser.add_argument('--steps',type=int,default=60)
    parser.add_argument('--duration',type=float,default=1e-4);parser.add_argument('--atol',type=float,default=1e-12)
    args=parser.parse_args()
    run(args.destination,json.loads(args.spec),epsilon=args.epsilon,steps=args.steps,
        duration=args.duration,absolute_tolerance=args.atol)
