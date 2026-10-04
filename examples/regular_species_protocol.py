"""Immutable species verification CLI with original-platform bitwise replay."""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
from plasma_bench.regular_species_transport import FrozenSpeciesTransport
import examples.regular_species_study as rates
import examples.regular_species_evolution as evolution
import examples.regular_species_controls as controls
import examples.species_shear_reference as shear
from .regular_species_records import provenance, sha, freeze


def protocol():
    cases = [('rate', dict(elements=n, degree=p)) for p in (3, 5) for n in (4, 8, 16, 32)]
    cases += [('rate', dict(elements=8, quadrature=q)) for q in (8, 16)]
    cases += [('rate', dict(elements=8, density=n)) for n in (1e16, 1e20)]
    cases += [('rate', dict(elements=8, contrast=0.)), ('rate', dict(elements=8, field_nr=32))]
    cases += [('evolution', dict(elements=n)) for n in (4, 8, 16)]
    cases += [('evolution', dict(elements=8, steps=n)) for n in (120, 240, 960)]
    cases += [('evolution', dict(elements=8, sign=-1)), ('evolution', dict(elements=8, species=1))]
    cases += [('control', dict(elements=n, points=8, field_nz=8)) for n in (2, 4)]
    return cases


def run(directory, selected=None):
    directory = Path(directory)
    identity = provenance()
    if directory.exists():
        if json.loads((directory/'provenance.json').read_text()) != identity:
            raise ValueError('Frozen source/runtime mismatch')
    else:
        directory.mkdir(parents=True)
        freeze(directory, identity)
    cases = protocol()
    for i in range(len(cases)) if selected is None else selected:
        kind, spec = cases[i]
        path = directory/f'case_{i:02d}.json'
        raw = path.with_suffix('.npz')
        if path.exists() or raw.exists():
            raise ValueError('Refusing to overwrite terminal or partial evidence')
        try:
            result, arrays = {'rate': rates, 'evolution': evolution, 'control': controls}[kind].run_case(spec)
            np.savez_compressed(raw, **arrays)
            result.update(arrays_file=raw.name, arrays_sha256=sha(raw))
        except Exception as exc:
            result = dict(spec=spec, status='failed', error=repr(exc), traceback=traceback.format_exc())
        result.update(kind=kind, source_sha256=identity['source_sha256'])
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(path.name, kind, spec, result['status'], flush=True)
        if result['status'] != 'complete':
            raise RuntimeError('Failed case preserved; no automatic retry')
    if provenance() != identity:
        raise RuntimeError('Source changed during study')


def audit(directory):
    directory = Path(directory)
    identity = json.loads((directory/'provenance.json').read_text())
    if provenance() != identity:
        raise ValueError('Original source/runtime required for bitwise replay')
    if (directory/'audit.json').exists():
        raise ValueError('Terminal audit already exists')
    records = []
    for i, (kind, spec) in enumerate(protocol()):
        path = directory/f'case_{i:02d}.json'
        record = json.loads(path.read_text())
        assert record['spec'] == spec and record['kind'] == kind and record['status'] == 'complete'
        assert record['source_sha256'] == identity['source_sha256']
        raw = directory/record['arrays_file']
        assert sha(raw) == record['arrays_sha256']
        arrays = dict(np.load(raw))
        replay, repeated = {'rate': rates, 'evolution': evolution, 'control': controls}[kind].run_case(spec)
        assert set(repeated) == set(arrays)
        for key, value in repeated.items():
            assert np.array_equal(value, arrays[key]), (path.name, key)
        for key, value in replay.items():
            if key != 'elapsed_seconds':
                assert value == record[key], (path.name, key)
        result = dict(file=path.name, kind=kind, spec=spec, arrays_and_metrics_bitwise=True)
        if kind in ('rate', 'evolution'):
            module = rates if kind == 'rate' else evolution
            result['order16'] = module.audit(spec, arrays, 16)
            result['order20'] = module.audit(spec, arrays, 20)
        if kind == 'evolution':
            _, _, _, model, _ = shear.initialize(spec)
            stepper = FrozenSpeciesTransport(model, arrays['current'], spec.get('species', 0))
            assert np.array_equal(stepper.step(arrays['penultimate_entropy'], float(arrays['dt'])), arrays['final_entropy'])
            result['final_step_bitwise'] = True
        records.append(result)
        (directory/'audit_progress.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], terminal=False, records=records), indent=2)+'\n')
        print(path.name, 'bitwise and independent audit complete', flush=True)
    assert provenance() == identity
    (directory/'audit.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], terminal=True, records=records), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--audit', action='store_true')
    parser.add_argument('--cases', type=int, nargs='*')
    args = parser.parse_args()
    audit(args.directory) if args.audit else run(args.directory, args.cases)
