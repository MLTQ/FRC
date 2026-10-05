"""Immutable physical-mode, boundary-control and full wall-RHS diagnostics."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.strong_wall_coordinates import StrongWallCoordinates
from plasma_bench.strong_wall_preconditioner import magnetic_left
from .strong_wall_reference import initialize
from .contact_mode_fields import reconstruct
from .contact_mode_projection import project
from .implicit_plasma_records import provenance, freeze, sha


def restricted(matrix, c, state, normal=True):
    nf = c.full.nfluid; keep = c.fluid_keep if normal else np.arange(nf)
    q = c.magnetic_basis; left = magnetic_left(c, state)
    return np.block([[matrix[np.ix_(keep, keep)], matrix[keep, nf:]@q],
        [left@matrix[nf:, keep], left@matrix[nf:, nf:]@q]])


def calculate(profile, epsilon, elements=1, points=8, field_nr=6, field_nz=8):
    spec = dict(kind='contact', elements=elements, points=points, field_nr=field_nr, field_nz=field_nz)
    _, _, u, system, b = initialize(spec); u[..., 1:4] = 0.
    c = StrongWallCoordinates(system, u, b)
    if profile == 'uniform': u = c.background()
    elif profile != 'contact': raise ValueError('Known stationary profile required')
    magnetic = b.packed_increment(system.geometry); maximum = 0.; wall = 0.
    def rate(value):
        nonlocal maximum, wall
        du, db = c.unpack(value); result = system.derivative(u+du, c.field(magnetic+db))
        maximum = max(maximum, result.residual)
        wall = max(wall, system.geometry.wall_residual(system.geometry.unpack(result.magnetic)[1]))
        return c.pack(result.fluid, result.magnetic)
    base = rate(np.zeros(c.size)); columns = []
    for k in range(c.size):
        direction = np.zeros(c.size); direction[k] = epsilon
        columns.append((rate(direction)-rate(-direction))/(2*epsilon))
    matrix = np.column_stack(columns)
    values, vectors = np.linalg.eig(matrix); index = int(values.real.argmax())
    vector = vectors[:, index]
    direction = vector.real if np.linalg.norm(vector.real) >= np.linalg.norm(vector.imag) else vector.imag
    direction = direction/np.max(abs(direction))
    mixed = (rate(epsilon*direction)-rate(-epsilon*direction))/(2*epsilon)
    predicted = matrix@direction
    result = dict(profile=profile, epsilon=epsilon, spec=spec, size=c.size,
        maximum_real=float(values.real[index]), leading_imaginary=float(values.imag[index]),
        spectral_radius=float(np.max(abs(values))), stationary_rate_max=float(np.max(abs(base))),
        mixed_direction_relative=float(np.linalg.norm(mixed-predicted)/np.linalg.norm(predicted)),
        eigenpair_backward=float(np.linalg.norm(matrix@vector-values[index]*vector)/(np.linalg.norm(matrix)*np.linalg.norm(vector))),
        max_induction_residual=maximum, max_rate_wall_residual=wall)
    arrays = dict(operator=matrix, eigenvalues=values, leading_vector=vector, base_rate=base,
        state=u, magnetic=magnetic, fluid_keep=c.fluid_keep, magnetic_basis=c.magnetic_basis,
        mixed_direction=direction, mixed_rate=mixed)
    return result, arrays


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.contact_wall_diagnosis', 'tests.test_strong_wall', 'tests.test_contact_mode_fields']
    identity = provenance(seeds); destination.mkdir(parents=True); freeze(destination, identity)
    physical = []; projections = []
    for order in (16, 20):
        result, fields = reconstruct(order); raw = destination/f'physical_{order}.npz'
        np.savez_compressed(raw, **fields); result.update(arrays_file=raw.name, arrays_sha256=sha(raw)); physical.append(result)
        for points, modes in ((32, 16), (64, 8), (64, 16), (64, 32)):
            projected = project(fields, points, modes); weight = projected['radial_norm'][:, None, None]
            value = projected['projected']; difference = projected['predicted']-result['eigenvalue']*value
            denominator = np.sum(weight*(result['eigenvalue']*value)**2, axis=(0, 2))
            raw = destination/f'projection_{order}_{points}_{modes}.npz'; np.savez_compressed(raw, **projected)
            projections.append(dict(order=order, points=points, modes=modes,
                captured_norm_fraction=(projected['captured_norms']/projected['original_norms']).tolist(),
                component_eigen_residual=np.sqrt(np.sum(weight*difference**2, axis=(0, 2))/denominator).tolist(),
                arrays_file=raw.name, arrays_sha256=sha(raw)))
        print('physical mode reconstructed', order, flush=True)
    original = Path('outputs/stationary_contact_linear_probe'); record = json.loads((original/'case_02.json').read_text())
    raw = original/record['arrays_file']
    if sha(raw) != record['arrays_sha256']: raise ValueError('Frozen original Jacobian required')
    original_matrix = np.load(raw)['operator']
    _, _, u, system, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    u[..., 1:4] = 0.; c = StrongWallCoordinates(system, u, b)
    nf = c.full.nfluid; fluid_only_keep = np.r_[c.fluid_keep, np.arange(nf, c.full.size)]
    controls = dict(normal_velocity_only=original_matrix[np.ix_(fluid_only_keep, fluid_only_keep)],
        magnetic_trace_only=restricted(original_matrix, c, u, False), both=restricted(original_matrix, c, u))
    raw = destination/'restricted_controls.npz'; np.savez_compressed(raw, **controls)
    control_metrics = {k: float(np.linalg.eigvals(a).real.max()) for k, a in controls.items()}
    cases = []
    for i, (profile, epsilon) in enumerate((('uniform', 1e-5), ('contact', 1e-4), ('contact', 1e-5))):
        start = time.monotonic(); result, arrays = calculate(profile, epsilon)
        if profile == 'contact':
            result['original_variational_matrix_relative'] = float(np.linalg.norm(arrays['operator']-controls['both'])/np.linalg.norm(controls['both']))
        raw_case = destination/f'case_{i:02}.npz'; np.savez_compressed(raw_case, **arrays)
        result.update(arrays_file=raw_case.name, arrays_sha256=sha(raw_case), elapsed_seconds=time.monotonic()-start)
        cases.append(result); (destination/f'case_{i:02}.json').write_text(json.dumps(result, indent=2)+'\n')
        print('wall Jacobian complete', profile, epsilon, result['maximum_real'], flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], physical_mode=physical, projections=projections,
        restricted_maximum_real=control_metrics, restricted_arrays_sha256=sha(raw), cases=cases,
        interpretation='Bounded-wall candidate. Original mode has incompatible physical wall traces; slow Bessel capture limits a matched continuum conclusion. Restricted matrices are diagnostic; separately formed full RHS Jacobians verify the actual candidate. Local current/enthalpy interfaces and refined full trajectories remain open.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    import sys
    run(sys.argv[1])
