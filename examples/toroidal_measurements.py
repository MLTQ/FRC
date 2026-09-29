"""Reproducible nonlinear toroidal evolution, conservative budgets and continuum errors."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import platform
import sys
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from scipy.special import j1, jn_zeros
from plasma_bench.toroidal_state import primitive, radius_squared
from .toroidal_cases import protocol
from .toroidal_field_errors import physical_errors
from .verification_records import write_record


def provenance():
    names = [str(p) for p in sorted(Path('plasma_bench').glob('toroidal_*.py'))]
    names += [f'plasma_bench/{name}.py' for name in ('mhd', 'mhd_mesh', 'mhd_state', 'mhd_flux', 'mhd_induction', 'mhd_reconstruction', 'mhd_positivity', 'magnetic_quadrature')]
    names += [f'examples/{name}.py' for name in ('toroidal_references', 'toroidal_cases', 'toroidal_measurements', 'toroidal_field_errors', 'verification_records')]
    hashes = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}
    return {'source_sha256': hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(), 'source_files': hashes,
            'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform()}


def field_errors(mesh, state, reference, gamma):
    q, exact = primitive(mesh, state, gamma), primitive(mesh, reference, gamma)
    volume, r2 = mesh.volumes, radius_squared(mesh)
    result = {}
    for name, index, weight in (('density', 0, volume), ('pressure', 3, volume),
                                ('rotation', 6, volume*r2), ('toroidal_field', 7, volume*r2)):
        diff = np.sum(weight*(q[..., index]-exact[..., index])**2)
        base = np.sum(weight*exact[..., index]**2)
        result[name+'_relative_L2'] = float(np.sqrt(diff/base)) if base > 0 else None
        result[name+'_absolute_L2'] = float(np.sqrt(diff/np.sum(weight)))
    angular_difference = state.conserved[..., 4]-reference.conserved[..., 4]
    result['angular_momentum_absolute_L2'] = float(np.sqrt(np.sum(volume*angular_difference**2)/np.sum(volume)))
    result['radial_velocity_rms_m_s'] = float(np.sqrt(np.sum(volume*q[..., 1]**2)/np.sum(volume)))
    return result


def measure(case, directory):
    source, start = provenance(), perf_counter()
    solver = case.build()
    mesh, initial = solver.mesh, solver.state.copy()
    reference = case.reference(mesh, case.end_time)
    history = {key: [] for key in ('time_s', 'total_energy_j', 'internal_energy_j', 'kinetic_energy_j', 'magnetic_energy_j',
                                  'toroidal_magnetic_energy_j', 'physical_ohmic_heat_j', 'mass_kg', 'angular_momentum_kg_m2_s',
                                  'toroidal_flux_wb', 'min_density_kg_m3', 'min_pressure_pa', 'max_face_divergence_t_m')}
    def record():
        data = solver.diagnostics()
        for key in history:
            history[key].append(data[key])
    record()
    initial_force, _ = solver.rhs(initial)
    while solver.time < case.end_time:
        remaining = case.end_time-solver.time
        if remaining <= 4*np.finfo(float).eps*case.end_time:
            solver.time = case.end_time
            break
        solver.step(min(remaining, solver.timestep_limit(), np.inf if case.dt_max is None else case.dt_max))
        record()
    final = solver.diagnostics()
    metrics = field_errors(mesh, solver.state, reference, case.gamma)
    physical6, physical10 = physical_errors(case, solver, 6), physical_errors(case, solver, 10)
    metrics['physical_reconstruction_relative_L2'] = physical10
    metrics['physical_error_quadrature_absolute_change'] = max(abs(physical6[k]-physical10[k]) for k in physical6 if physical6[k] is not None)
    if case.kind == 'decay':
        k = jn_zeros(0, 1)[0]/mesh.radius
        continuum_energy0 = np.pi*(2*mesh.half_length)/MU0*case.amplitude**2*(mesh.radius**2/2)*j1(k*mesh.radius)**2
        exact_heat = continuum_energy0*(-np.expm1(-2*case.diffusivity*k*k*case.end_time))
        metrics['continuum_heat_j'] = float(exact_heat)
        metrics['physical_heat_relative_error'] = float((final['physical_ohmic_heat_j']-exact_heat)/exact_heat)
        metrics['axis_g_relative_error'] = float((solver.state.conserved[0, 0, 5]-reference.conserved[0, 0, 5])/reference.conserved[0, 0, 5])
    arrays = {'initial_conserved': initial.conserved, 'initial_psi': initial.psi,
              'final_conserved': solver.state.conserved, 'final_psi': solver.state.psi,
              'reference_conserved': reference.conserved, 'reference_psi': reference.psi,
              'initial_rhs_conserved': initial_force.conserved, 'initial_rhs_psi': initial_force.psi,
              'r_faces': mesh.rf, 'z_faces': mesh.zf, 'volumes': mesh.volumes}
    path = directory/(case.name+'.npz')
    np.savez_compressed(path, **arrays)
    result = {'completed': True, 'name': case.name, 'configuration': case.configuration(), 'provenance': source,
              'elapsed_seconds': perf_counter()-start, 'initial': solver.initial_integrals, 'final': final,
              'errors': metrics, 'history_columns': history,
              'fields': {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'units': 'SI; six conserved variables documented in toroidal_state; psi Wb/(2pi), coordinates m'}}
    if provenance() != source:
        raise RuntimeError('Toroidal source changed during measurement')
    write_record(directory/(case.name+'.json'), result)
    print('Completed', case.name, final['steps'], 'steps', flush=True)
    return result


def run(directory='outputs/toroidal_verification'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    cases = [measure(case, directory) for case in protocol()]
    if provenance() != source:
        raise RuntimeError('Source changed across toroidal study')
    write_record(directory/'measurements.json', {'completed': True, 'provenance': source, 'cases': cases})


if __name__ == '__main__':
    run()
