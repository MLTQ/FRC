"""Run and assemble the complete conservative warm Hall verification protocol."""
import argparse
from pathlib import Path
from time import perf_counter
import hashlib
import json
from .warm_hall_cases import protocol
from .warm_hall_records import provenance, save
from .warm_hall_spatial import pressure_spatial, trace_diffusion
from .warm_pressure_evolution import measure as pressure_evolution
from .warm_hall_evolution import measure as coupled_evolution
from .verification_records import write_record


def group(case):
    if case['kind'].startswith('pressure'):
        return 'pressure'
    if case['kind'] == 'trace_diffusion':
        return 'trace'
    if case['kind'] == 'nonlinear':
        return 'nonlinear'
    return 'wave_slow' if case['reference']['branch'] == 0 else 'wave_other'


def assemble(directory):
    source, records = provenance(), []
    for case in protocol():
        record = json.loads((directory/(case['name']+'.json')).read_text())
        assert record['completed'] and record['configuration'] == case and record['provenance'] == source
        assert hashlib.sha256((directory/record['fields']['file']).read_bytes()).hexdigest() == record['fields']['sha256']
        records.append(record)
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, cases=records))
    print('Assembled', len(records), source['source_sha256'], flush=True)


def run(directory, selected=None):
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    for case in protocol():
        if selected is not None and group(case) != selected:
            continue
        start = perf_counter()
        fn = {'pressure_spatial': pressure_spatial, 'trace_diffusion': trace_diffusion,
              'pressure_evolve': pressure_evolution, 'wave': coupled_evolution, 'nonlinear': coupled_evolution}[case['kind']]
        record, arrays = fn(case)
        save(directory, case, source, record | dict(elapsed_seconds=perf_counter()-start), arrays)
        print('Completed', case['name'], record['metrics'], flush=True)
    if selected is None:
        assemble(directory)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--group', choices=['pressure', 'trace', 'wave_slow', 'wave_other', 'nonlinear'])
    parser.add_argument('--assemble', action='store_true')
    parser.add_argument('--directory', default='outputs/warm_hall_verification')
    args = parser.parse_args()
    if args.assemble:
        assemble(Path(args.directory))
    else:
        run(Path(args.directory), args.group)
