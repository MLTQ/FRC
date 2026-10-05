"""Complete immutable replay of signed nonlinear trajectory checkpoints."""
import argparse
import json
from pathlib import Path
from .matched_mode_evolution import run as ordinary_run
from .guarded_mode_evolution import run as guarded_run
from .wall_diagnostic_replay import verify


def run(original,repeated,destination,*,guarded=True):
    original=Path(original);record=json.loads((original/'summary.json').read_text())
    if not record['terminal'] or not record['success']:raise ValueError('Complete passing original trajectory required')
    function=guarded_run if guarded else ordinary_run
    repeated_record=function(repeated,record['spec'],epsilon=record['epsilon'],steps=record['steps'],
        duration=record['duration'],absolute_tolerance=record['absolute_tolerance'])
    if not repeated_record['success']:raise RuntimeError('Replay trajectory rejected; no successful replay claimed')
    verify(original,repeated,destination)
    print('complete trajectory bitwise verified',original,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('original');parser.add_argument('repeated');parser.add_argument('destination')
    parser.add_argument('--unguarded',action='store_true');args=parser.parse_args()
    run(args.original,args.repeated,args.destination,guarded=not args.unguarded)
