"""Explicit finer-time/mesh and slow-wave sensitivity controls, separate protocol."""
import argparse
from copy import deepcopy
from pathlib import Path
import hashlib
import json
from time import perf_counter
from .warm_hall_cases import protocol as base_protocol
from .warm_hall_records import provenance, save
from .warm_pressure_evolution import measure as pressure_measure
from .warm_hall_evolution import measure as wave_measure
from .verification_records import write_record


def protocol():
    base = {c['name']: c for c in base_protocol()}; cases = []
    for steps in (160, 320, 640):
        c = deepcopy(base['pressure_evolve_64'])
        c.update(name=f'pressure_fine_time_{steps}', steps=steps); cases.append(c)
    c = deepcopy(base['pressure_evolve_64'])
    c.update(name='pressure_fine_mesh_128', nr=128, nz=256, steps=320); cases.append(c)
    c = deepcopy(base['wave_b0_16'])
    c.update(name='slow_wave_half_cfl', fluid_cfl=.125); cases.append(c)
    c = deepcopy(base['wave_b0_16'])
    c.update(name='slow_wave_amplitude_01'); c['reference']['amplitude'] = .1; cases.append(c)
    return cases


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def run(group=None, directory='outputs/warm_hall_controls'):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source, control = provenance(), source_hash()
    for case in protocol():
        label = 'pressure' if case['kind'] == 'pressure_evolve' else 'slow'
        if group is not None and group != label:
            continue
        start = perf_counter()
        record, arrays = (pressure_measure if label == 'pressure' else wave_measure)(case)
        assert source_hash() == control
        save(directory, case, source, record | dict(control_source_sha256=control, elapsed_seconds=perf_counter()-start), arrays)
        print('Completed', case['name'], record['metrics'], flush=True)
    if group is None:
        assemble(directory)


def assemble(directory='outputs/warm_hall_controls'):
    directory = Path(directory); records = []; source = provenance()
    for case in protocol():
        r = json.loads((directory/(case['name']+'.json')).read_text())
        assert r['completed'] and r['configuration'] == case and r['provenance'] == source
        assert r['control_source_sha256'] == source_hash()
        assert hashlib.sha256((directory/r['fields']['file']).read_bytes()).hexdigest() == r['fields']['sha256']
        records.append(r)
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, control_source_sha256=source_hash(), cases=records))
    print('Assembled controls', len(records), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--group', choices=['pressure', 'slow'])
    parser.add_argument('--assemble', action='store_true')
    args = parser.parse_args()
    assemble() if args.assemble else run(args.group)
