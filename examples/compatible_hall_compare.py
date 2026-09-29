"""Compare archived Hall trajectories without changing their evolved states."""
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.compatible_magnetic import energy_density
from plasma_bench.toroidal_state import radius_squared
from .compatible_hall_records import relative


def state_difference(mesh, actual, reference, initial):
    """Normalize differences by the reference change, never the uniform bias."""
    u, psi = actual
    ref, pref = reference
    start, pstart = initial
    weight, s = mesh.volumes, radius_squared(mesh)
    def magnetic_norm(p, g):
        return float(np.sum(weight*(energy_density(mesh, p)+s*g*g/(2*MU0))))
    numerator = magnetic_norm(psi-pref, u[..., 5]-ref[..., 5])
    denominator = magnetic_norm(pref-pstart, ref[..., 5]-start[..., 5])
    result = dict(magnetic_relative_to_reference_change=np.sqrt(numerator/denominator) if denominator else None)
    for label, column in (('density', 0), ('energy', 3)):
        result[label+'_relative_to_reference_change'] = relative(u[..., column]-ref[..., column], ref[..., column]-start[..., column], weight)
    # Convert angular momentum to azimuthal momentum using the exact cell r² moment.
    scale = np.stack(np.broadcast_arrays(np.ones_like(s), np.ones_like(s), 1/np.sqrt(s)), axis=-1)
    columns = [1, 2, 4]
    result['momentum_relative_to_reference_change'] = relative((u[..., columns]-ref[..., columns])*scale,
                                                             (ref[..., columns]-start[..., columns])*scale, weight[..., None])
    return result


def restrict(fine, coarse_mesh):
    """Exact annular-volume restriction of cells; coincident corner flux."""
    u, psi, volume = fine['final_conserved'], fine['final_psi'], fine['volumes']
    nr, nz = coarse_mesh.nr, coarse_mesh.nz
    rr, rz = u.shape[0]//nr, u.shape[1]//nz
    if u.shape[:2] != (nr*rr, nz*rz):
        raise ValueError('Nested integer meshes required')
    restricted = (u*volume[..., None]).reshape(nr, rr, nz, rz, 6).sum(axis=(1, 3))/coarse_mesh.volumes[..., None]
    return restricted, psi[::rr, ::rz]


def comparisons(records, archives, meshes):
    time = {}
    for family, counts, finest in (('wave', (8, 16, 32), 256), ('nonlinear', (4, 8, 16), 128)):
        ref = archives[f'{family}_time_{finest}']
        for count in counts:
            name = f'{family}_time_{count}'
            a = archives[name]
            time[name] = state_difference(meshes[name], (a['final_conserved'], a['final_psi']),
                                          (ref['final_conserved'], ref['final_psi']),
                                          (ref['initial_conserved'], ref['initial_psi']))
            time[name]['actual_steps'] = records[name]['final']['steps']
    spatial = {}
    ref = archives['nonlinear_n64']
    for nr in (8, 16, 32):
        name = f'nonlinear_n{nr}'
        a, mesh = archives[name], meshes[name]
        initial_ref = dict(final_conserved=ref['initial_conserved'], final_psi=ref['initial_psi'], volumes=ref['volumes'])
        spatial[name] = state_difference(mesh, (a['final_conserved'], a['final_psi']), restrict(ref, mesh), restrict(initial_ref, mesh))
    return dict(fixed_grid_time=time, nonlinear_nested_mesh=spatial,
                normalization='Difference divided by fine-reference change from its own initial state; magnetic quadratic compatible norm, fluid volume L2. Mesh comparison restricts conservative cells and corner psi, not a continuum physical-field norm.')


def wave_modes(records, archives, meshes):
    """Resolve fundamental amplitude/phase using independent exact cell templates."""
    from .compressible_hall_reference import CompressibleHallWave
    from .compatible_hall_references import wave_state
    results = {}
    for name, record in records.items():
        case = record['configuration']
        if case['kind'] != 'wave' or not case['reference']['hall']:
            continue
        mesh = meshes[name]
        ref = CompressibleHallWave(**case['reference'])
        c = wave_state(mesh, ref).conserved[..., 5]
        s = wave_state(mesh, ref, np.pi/(2*abs(ref.omega))).conserved[..., 5]
        basis = np.stack((c, s), axis=-1)
        weight = mesh.volumes*radius_squared(mesh)
        gram = np.einsum('ij,ijk,ijl->kl', weight, basis, basis)
        actual = archives[name]['final_conserved'][..., 5]
        coefficients = np.linalg.solve(gram, np.einsum('ij,ijk,ij->k', weight, basis, actual))
        phase = np.arctan2(coefficients[1], coefficients[0])
        error = np.arctan2(np.sin(phase-case['phase']), np.cos(phase-case['phase']))
        fit = np.einsum('ijk,k->ij', basis, coefficients)
        results[name] = dict(fundamental_amplitude=float(np.linalg.norm(coefficients)), phase_error_radians=float(error),
                            relative_phase_error=float(error/case['phase']),
                            remaining_g_relative_L2=relative(actual-fit, actual, weight))
    return results
