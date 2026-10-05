"""Frozen complete moving wave/rate/spectrum protocol and independent replay."""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
import examples.regular_plasma_rates as rates
import examples.regular_plasma_evolution as evolution
import examples.regular_plasma_linear_spectrum as spectrum
from plasma_bench.regular_plasma import RegularPlasmaRate
from plasma_bench.regular_plasma_step import field_from_packed, rk4
from .regular_plasma_reference import initialize
from .regular_plasma_norms import measure
from .regular_plasma_records import provenance, sha, freeze


def protocol():
    cases = [('rate', dict(kind=k, elements=n)) for k in ('contact', 'wave') for n in (4, 8, 16)]
    cases += [('rate', dict(kind='wave', elements=8, amplitude=a)) for a in (1., .01)]
    cases += [('rate', dict(kind='wave', elements=8, branch=b)) for b in (0, 1)]
    cases += [('rate', dict(kind='wave', elements=8, quadrature=q)) for q in (8, 16)]
    cases += [('evolution', dict(kind='wave', elements=n, amplitude=1., steps=120)) for n in (4, 8, 16)]
    cases += [('evolution', dict(kind='wave', elements=4, amplitude=1., steps=n)) for n in (60, 240, 480)]
    cases += [('evolution', dict(kind='wave', elements=8, amplitude=.1, steps=120))]
    cases += [('spectrum', dict(elements=1, points=8, field_nr=6, field_nz=8, epsilon=e)) for e in (1e-4, 1e-5)]
    cases += [('rate', dict(kind='wave', elements=n, branch=b)) for b in (0, 1) for n in (16, 32)]
    return cases


def run(directory, selected=None):
    directory = Path(directory); identity = provenance()
    if directory.exists():
        if json.loads((directory/'provenance.json').read_text()) != identity:
            raise ValueError('Frozen source/runtime mismatch')
    else:
        directory.mkdir(parents=True); freeze(directory, identity)
    cases = protocol()
    for i in range(len(cases)) if selected is None else selected:
        kind, spec = cases[i]; path = directory/f'case_{i:02d}.json'; raw = path.with_suffix('.npz')
        if path.exists() or raw.exists():
            raise ValueError('Refusing to overwrite terminal or partial evidence')
        try:
            result, arrays = {'rate': rates, 'evolution': evolution, 'spectrum': spectrum}[kind].run_case(spec)
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
    directory = Path(directory); identity = json.loads((directory/'provenance.json').read_text())
    if provenance() != identity:
        raise ValueError('Original source/runtime required for bitwise replay')
    if (directory/'audit.json').exists():
        raise ValueError('Terminal audit already exists')
    records = []
    for i, (kind, spec) in enumerate(protocol()):
        path = directory/f'case_{i:02d}.json'; record = json.loads(path.read_text())
        assert record['spec'] == spec and record['kind'] == kind and record['status'] == 'complete'
        assert record['source_sha256'] == identity['source_sha256']
        raw = directory/record['arrays_file']; assert sha(raw) == record['arrays_sha256']
        arrays = dict(np.load(raw))
        replay, repeated = {'rate': rates, 'evolution': evolution, 'spectrum': spectrum}[kind].run_case(spec)
        assert set(repeated) == set(arrays)
        for key, value in repeated.items():
            assert np.array_equal(value, arrays[key]), (path.name, key)
        for key, value in replay.items():
            if key != 'elapsed_seconds':
                assert value == record[key], (path.name, key)
        result = dict(file=path.name, kind=kind, spec=spec, arrays_and_metrics_bitwise=True)
        if kind != 'spectrum':
            ref, d, u, system, field = initialize(spec)
            for order in (16, 20):
                if kind == 'rate':
                    rate = RegularPlasmaRate(arrays['fluid_rate'], arrays['magnetic_rate'],
                        arrays['current_rate'], arrays['electric'], record['residual'], record['diagnostics'])
                    result[f'order{order}'] = measure(spec, ref, d, u, system, field, rate, order)
                else:
                    final_field = field_from_packed(system.geometry, arrays['final_magnetic'], field.axial_bias)
                    result[f'order{order}'] = evolution.field_norms(spec, ref, d, arrays['final_state'], system,
                                                                  final_field, record['duration'], order)
            if kind == 'evolution':
                previous = field_from_packed(system.geometry, arrays['penultimate_magnetic'], field.axial_bias)
                final, b, _ = rk4(system, arrays['penultimate_state'], previous, float(arrays['dt']))
                assert np.array_equal(final, arrays['final_state'])
                assert np.array_equal(b.packed_increment(system.geometry), arrays['final_magnetic'])
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
