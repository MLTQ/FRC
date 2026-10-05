"""Frozen stationary-contact Jacobian diagnostic, separate from time stepping."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.refined_regular_plasma import RefinedRegularPlasma
from plasma_bench.plasma_coordinates import PlasmaCoordinates
from .regular_plasma_reference import initialize
from .implicit_plasma_records import provenance, freeze, sha


def calculate(profile, epsilon):
    ref, d, u, old, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    u[..., 1:4] = 0.
    system = RefinedRegularPlasma(old.mapping, ref.ion_mass, **old.solver_options)
    c = PlasmaCoordinates(system, u, b)
    if profile == 'uniform': u = c.background()
    elif profile != 'contact': raise ValueError('Known stationary profile required')
    magnetic = b.packed_increment(system.geometry)
    max_residual = 0.; max_refinement = 0
    def rate(value):
        nonlocal max_residual, max_refinement
        du, db = c.unpack(value); result = system.derivative(u+du, c.field(magnetic+db))
        max_residual = max(max_residual, result.residual)
        max_refinement = max(max_refinement, result.diagnostics['induction_residuals'].get('outer_refinements', 0))
        return c.pack(result.fluid, result.magnetic)
    base = rate(np.zeros(c.size)); columns = []
    for k in range(c.size):
        direction = np.zeros(c.size); direction[k] = epsilon
        columns.append((rate(direction)-rate(-direction))/(2*epsilon))
    matrix = np.column_stack(columns)
    values, vectors = np.linalg.eig(matrix); index = int(values.real.argmax())
    direction = vectors[:, index]
    direction = direction.real if np.linalg.norm(direction.real) >= np.linalg.norm(direction.imag) else direction.imag
    direction /= np.max(abs(direction))
    mixed = (rate(epsilon*direction)-rate(-epsilon*direction))/(2*epsilon)
    predicted = matrix@direction
    result = dict(profile=profile, epsilon=epsilon, size=c.size, bias=b.axial_bias,
        maximum_real=float(values.real[index]), leading_imaginary=float(values.imag[index]),
        spectral_radius=float(np.max(abs(values))), stationary_rate_max=float(np.max(abs(base))),
        eigenpair_backward=float(np.linalg.norm(matrix@vectors[:, index]-values[index]*vectors[:, index])/
            (np.linalg.norm(matrix)*np.linalg.norm(vectors[:, index]))),
        mixed_direction_relative_error=float(np.linalg.norm(mixed-predicted)/np.linalg.norm(predicted)),
        max_induction_residual=max_residual, max_refinements=max_refinement)
    return result, dict(state=u, magnetic=magnetic, coordinates=c.scale, operator=matrix,
        eigenvalues=values, leading_vector=vectors[:, index], base_rate=base,
        mixed_direction=direction, mixed_rate=mixed)


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.contact_linear_probe']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity); records = []
    for i, (profile, epsilon) in enumerate([('uniform', 1e-5), ('contact', 1e-4), ('contact', 1e-5)]):
        start = time.monotonic(); result, arrays = calculate(profile, epsilon)
        raw = destination/f'case_{i:02}.npz'; np.savez_compressed(raw, **arrays)
        result.update(arrays_file=raw.name, arrays_sha256=sha(raw), source_sha256=identity['source_sha256'], elapsed_seconds=time.monotonic()-start)
        (destination/f'case_{i:02}.json').write_text(json.dumps(result, indent=2)+'\n')
        records.append(result); print('complete', profile, epsilon, result['maximum_real'], flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Probe source changed')
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], terminal=True, cases=records), indent=2)+'\n')


if __name__ == '__main__': run('outputs/stationary_contact_linear_probe')
