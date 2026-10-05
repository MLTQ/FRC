"""Replay frozen spatial probe and compare separate continuum-mode diagnostics."""
import json
from pathlib import Path
import numpy as np
from plasma_bench.refined_regular_plasma import RefinedRegularPlasma
from plasma_bench.plasma_coordinates import PlasmaCoordinates
from .regular_plasma_reference import initialize
from .regular_plasma_linear import StationaryLinearPlasma
from .contact_linear_probe import calculate
from .stationary_contact_continuum import operator, uniform_wave_errors
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.stationary_contact_audit']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity)
    original = Path('outputs/stationary_contact_linear_probe')
    records = json.loads((original/'summary.json').read_text())['cases']; replayed = []
    for record in records:
        raw = original/record['arrays_file']
        if sha(raw) != record['arrays_sha256']: raise ValueError('Original raw hash required')
        repeated, arrays = calculate(record['profile'], record['epsilon'])
        expected = dict(np.load(raw))
        if any(not np.array_equal(value, expected[key]) for key, value in arrays.items()):
            raise RuntimeError('Spatial probe array replay failed')
        if any(record[key] != value for key, value in repeated.items()): raise RuntimeError('Spatial probe metric replay failed')
        replayed.append(dict(file=record['arrays_file'], arrays_sha256=sha(raw), bitwise_replayed=True))
    ref, _, u, old, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    u[..., 1:4] = 0.; system = RefinedRegularPlasma(old.mapping, ref.ion_mass, **old.solver_options)
    coordinates = PlasmaCoordinates(system, u, b); reference = StationaryLinearPlasma(system, coordinates.background(), b.axial_bias)
    uniform = np.load(original/'case_00.npz')['operator']; nf = coordinates.nfluid
    reference_errors = []
    for direction in np.random.default_rng(208).normal(size=(3, coordinates.size)):
        fluid, magnetic, _ = reference.derivative(*coordinates.unpack(direction))
        expected = coordinates.pack(fluid, magnetic)
        reference_errors.append(float(np.linalg.norm(uniform@direction-expected)/np.linalg.norm(expected)))
    full = np.load(original/'case_02.npz')['operator']
    scalars = np.flatnonzero(np.isin(coordinates.fluid_indices%6, [0, 4, 5]))
    entropy = np.flatnonzero(np.isin(coordinates.fluid_indices%6, [4, 5])); magnetic = np.arange(nf, coordinates.size)
    ablated = full.copy(); ablated[np.ix_(magnetic, scalars)] = 0.; ablated[np.ix_(entropy, magnetic)] = 0.
    fixed = np.r_[entropy, magnetic]
    components = {}
    for name, matrix in [('full', full), ('fluid_only', full[:nf, :nf]), ('magnetic_only', full[nf:, nf:]),
                         ('species_pair_removed', ablated), ('fixed_bulk_species', full[np.ix_(fixed, fixed)])]:
        values = np.linalg.eigvals(matrix)
        components[name] = dict(maximum_real=float(values.real.max()), spectral_radius=float(np.max(abs(values))))
    raw = destination/'component_controls.npz'; np.savez_compressed(raw, species_pair_removed=ablated, fixed_bulk_indices=fixed)
    continuum = []
    for i, (points, radial) in enumerate([(8, 1), (16, 1), (32, 1), (64, 1), (128, 1), (64, 4), (64, 8), (64, 16)]):
        matrix, arrays = operator(points, radial)
        values, vectors = np.linalg.eig(matrix); index = int(values.real.argmax())
        vector = vectors[:, index]; physical = (vector*arrays['scale']).reshape(9, points)
        divergence = arrays['q']*physical[6]+arrays['derivative']@physical[8]
        norm = np.linalg.norm(arrays['q']*physical[6])+np.linalg.norm(arrays['derivative']@physical[8])
        raw_case = destination/f'continuum_{i:02}.npz'
        np.savez_compressed(raw_case, **arrays, operator=matrix, eigenvalues=values, leading_vector=vector)
        continuum.append(dict(points=points, radial_mode=radial, q=float(arrays['q']),
            maximum_real=float(values.real[index]), leading_imaginary=float(values.imag[index]),
            divergence_relative=float(np.linalg.norm(divergence)/norm),
            eigenpair_backward=float(np.linalg.norm(matrix@vector-values[index]*vector)/
                (np.linalg.norm(matrix)*np.linalg.norm(vector))),
            arrays_file=raw_case.name, arrays_sha256=sha(raw_case)))
        print('continuum diagnostic', points, radial, values[index], flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Audit source changed')
    result = dict(source_sha256=identity['source_sha256'], terminal=True, probe_replay=replayed,
        uniform_independent_mixed_errors=reference_errors, component_controls=components,
        component_arrays_sha256=sha(raw), continuum_uniform_wave_errors=uniform_wave_errors(32), continuum=continuum,
        interpretation='Unmatched radial spaces: these are diagnostics, not a mode-by-mode proof that a production eigenvalue is unphysical. Uniform null-mode relative mixed error in the original probe is ill-conditioned; use the independent nonzero mixed-direction controls here.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__': run('outputs/stationary_contact_audit')
