"""Independent density reconstruction of an immutable accepted-prefix snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import shutil
import numpy as np
from scipy.special import j0
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def fingerprint(value):
    value = np.ascontiguousarray(value)
    header = json.dumps(dict(shape=value.shape, dtype=value.dtype.str), sort_keys=True).encode()
    return hashlib.sha256(header + value.tobytes()).hexdigest()


def calculate(case):
    directory = Path(case['directory']); record = case['summary_snapshot']
    steps = case['prefix_steps']; spec = record['spec']
    if type(steps) is not int or not 1 <= steps <= record['completed_steps']:
        raise ValueError('An existing nonempty accepted prefix is required')
    if set(spec) != {'elements', 'points', 'field_nr', 'field_nz'}:
        raise ValueError('Only the original four-count contact fixture is supported')
    source = json.loads((directory/'provenance.json').read_text())
    if source['source_sha256'] != record['source_sha256']:
        raise ValueError('Unchanged physical source identity required')
    if any(sha(p) != h for p,h in source['source_files'].items()):
        raise ValueError('Original physical calculation source required')
    inputs = {'initial.npz': record['initial_sha256']}
    selected = record['records'][:steps]
    if [r['step'] for r in selected] != list(range(1, steps+1)):
        raise ValueError('Every consecutive accepted endpoint is required')
    inputs.update({r['arrays_file']: r['arrays_sha256'] for r in selected})
    if any(sha(directory/p) != h for p,h in inputs.items()):
        raise ValueError('Saved initial and accepted state/stage bytes changed')
    mode = ContinuumContactMode(); ref = mode.ref
    if mode.input_sha256 != record['input_sha256'] or mode.growth != record['continuum_growth']:
        raise ValueError('The original independent continuum reference is required')
    fluid = RegularCoordinateFluid(ref.radius, ref.length, spec['elements'], spec['points'],
        ref.ion_mass+ref.electron_mass, degree=5, gammas=(ref.gamma_e, ref.gamma_i))
    mesh = MHDMesh(ref.radius, ref.length/2, spec['field_nr'], spec['field_nz'])
    if spec['field_nr'] % spec['elements']:
        raise ValueError('Quadrature cells must subdivide every fluid element')
    with np.load(directory/'initial.npz') as raw:
        initial = raw['state'][...,0].copy()
    if initial.shape != fluid.shape[:-1]: raise ValueError('Original density shape required')
    rows = []; arrays = {}; times = np.array([0., *[r['time'] for r in selected]])
    names = ['initial.npz', *[r['arrays_file'] for r in selected]]
    for order in (16,20):
        sample = TransportSampler(fluid, mesh, order)
        base = ref.fields(sample.Z)['number_density']
        direction = j0(mode.q*sample.R)*mode.amplitudes(sample.Z)[...,0]
        intended = record['epsilon']*direction
        initial_sample = sample.scalar(initial)
        initial_error = initial_sample-base-intended
        norm = lambda value: float(np.sqrt(np.sum(sample.weight*value*value)))
        signal_norm = norm(intended)
        if not np.isfinite(signal_norm) or signal_norm <= 0:
            raise ValueError('A finite nonzero reference density signal is required')
        metrics = []; hashes = []
        for filename, when in zip(names, times):
            with np.load(directory/filename) as raw:
                density = raw['state'][...,0].copy()
            if density.shape != initial.shape or not np.isfinite(density).all():
                raise ValueError('Finite unchanged-shape density required')
            # Subtract stored coefficients before independent linear reconstruction.
            delta = sample.scalar(density-initial)
            growth = float(np.exp(mode.growth*when))
            error = delta-np.expm1(mode.growth*when)*intended+initial_error
            sampled = initial_sample+delta
            metrics.append([norm(error)/(growth*signal_norm), norm(delta)/signal_norm,
                float(np.max(np.abs((density-initial)/initial))), float(density.min()),
                float(sampled.min()), float(np.max(np.abs(delta/initial_sample)))])
            hashes.append(dict(density_change=fingerprint(delta), reference_error=fingerprint(error)))
        arrays[f'q{order}_metrics'] = np.array(metrics)
        rows.append(dict(order=order, initial_signal_L2=signal_norm,
            initial_representation_error=norm(initial_error)/signal_norm,
            weight_sha256=fingerprint(sample.weight), reference_sha256=fingerprint(direction),
            physical_fingerprints=hashes))
    arrays['time'] = times; arrays['expected_growth'] = np.exp(mode.growth*times)
    if any(sha(directory/p) != h for p,h in inputs.items()):
        raise ValueError('Accepted prefix changed during observation')
    return dict(directory=str(directory), prefix_steps=steps, spec=spec,
        physical_source_sha256=record['source_sha256'], epsilon=record['epsilon'],
        original_terminal_at_snapshot=record['terminal'], original_success_at_snapshot=record['success'],
        input_raw_sha256=inputs, reference_input_sha256=mode.input_sha256, quadrature=rows), arrays


def run(manifest_path, destination):
    manifest_path = Path(manifest_path); destination = Path(destination)
    if destination.exists(): raise ValueError('Fresh density-observation destination required')
    manifest_hash = sha(manifest_path); manifest = json.loads(manifest_path.read_text())
    identity = provenance(['examples.trajectory_density_prefix'])
    destination.mkdir(parents=True); freeze(destination, identity); cases = []
    for case in manifest['cases']:
        record, arrays = calculate(case); path = destination/f'case_{len(cases):02}.npz'
        np.savez_compressed(path, **arrays)
        record.update(arrays_file=path.name, arrays_sha256=sha(path)); cases.append(record)
        print('accepted density prefix', record['directory'], record['prefix_steps'], flush=True)
    if sha(manifest_path) != manifest_hash or provenance(['examples.trajectory_density_prefix']) != identity:
        raise ValueError('Manifest or observation source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], manifest=str(manifest_path),
        manifest_sha256=manifest_hash, cases=cases,
        columns=['relative_L2_error_to_growing_linear_density_signal',
            'density_change_L2_over_initial_reference_signal', 'nodal_max_fractional_density_change',
            'nodal_density_min', 'sampled_density_min', 'sampled_max_fractional_density_change'],
        scope='Single-sign accepted prefixes, including any live trajectory. Density only, using the independent linear continuum target. Neither a full-duration pass, a centered linear response, a nonlinear exact solution nor causal identification. No production RHS or time integrator is called.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (destination/'resource_observation.json').write_text(json.dumps(dict(
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        free_bytes_after=shutil.disk_usage('.').free,
        scope='Own-process macOS cumulative peak and final disk observation; no continuous audit resource guard'),indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('manifest'); parser.add_argument('destination')
    args = parser.parse_args(); run(args.manifest,args.destination)
