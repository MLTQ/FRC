"""Reproducible actual-cylinder field, dispersion, dissipation and energy measurements."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import platform
import sys
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from .hall_wave_cases import protocol
from .hall_time_reference import exact_spatial_evolution
from .verification_records import write_record


def provenance():
    names = ['plasma_bench/hall_parameters.py', 'plasma_bench/hall_modes.py', 'plasma_bench/hall_cylinder.py',
             'examples/hall_wave_reference.py', 'examples/hall_wave_cases.py', 'examples/hall_wave_measurements.py',
             'examples/hall_time_reference.py', 'examples/verification_records.py']
    hashes = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}
    return {'source_files': hashes, 'source_sha256': hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
            'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform()}


def physical_checks(solver, reference, state, time, order):
    r, z, weights = solver.quadrature_points(order)
    fields, errors, energies = {}, {}, {}
    for quantity in ('magnetic', 'velocity', 'current'):
        measured = solver.sample(state, r, z, quantity)
        exact = reference.fields(r, z, time, quantity)
        fields[quantity] = measured
        denominator = np.sum(weights*np.sum(exact**2, axis=-1))
        errors[quantity+'_relative_volume_L2'] = float(np.sqrt(np.sum(weights*np.sum((measured-exact)**2, axis=-1))/denominator))
        energies[quantity] = float(np.sum(weights*np.sum(measured**2, axis=-1)))
    physical = {'magnetic': energies['magnetic']/(2*MU0),
                'bulk_kinetic': energies['velocity']*solver.parameters.density/2,
                'current_kinetic': energies['current']*MU0*solver.parameters.current_length_squared/2}
    heat = float(np.sum(weights*reference.heat_density(r, z, time)))
    return errors, physical, heat, {'r': r, 'z': z, 'volume_weights': weights, **fields}


def measure(case, directory):
    start = perf_counter()
    source = provenance()
    solver, reference, state = case.build()
    initial = state.copy()
    end = case.phase_duration/reference.omega.real
    dt = end/case.steps
    initial_energy = sum(solver.energies(state).values())
    history, heat = [], 0.
    mode_index = case.axial_mode
    template = initial[:, mode_index, :2]
    norm = np.vdot(template, template).real
    for step_number in range(1, case.steps+1):
        update = solver.step(state, dt)
        state, heat = update.state, heat+update.heat
        amplitude = np.vdot(template, state[:, mode_index, :2])/norm
        energy = solver.energies(state)
        history.append({'time_s': step_number*dt, 'energy_j': energy, 'heat_j': heat,
                        'relative_energy_balance': (sum(energy.values())+heat-initial_energy)/initial_energy,
                        'equation_relative_residual': update.equation_relative_residual,
                        'amplitude_real': float(amplitude.real), 'amplitude_imag': float(amplitude.imag)})
    errors, physical, exact_heat, fields = physical_checks(solver, reference, state, end, 8)
    errors12, physical12, exact_heat12, _ = physical_checks(solver, reference, state, end, 12)
    initial_errors, _, _, _ = physical_checks(solver, reference, initial, 0., 8)
    phase_error = -np.angle(amplitude*np.exp(1j*reference.omega.real*end))/case.phase_duration
    damping = -np.log(abs(amplitude))/end
    exact_spatial = exact_spatial_evolution(solver, initial, end)
    temporal_error = float(np.sqrt(sum(solver.modes.energies(state-exact_spatial, solver.volume_factor).values())/sum(solver.energies(exact_spatial).values())))
    arrays = {'initial_modes': initial, 'final_modes': state, 'exact_spatial_final_modes': exact_spatial, 'radial_nodes': solver.mesh.r,
              'axial_nodes': solver.mesh.z, 'radial_eigenvalues': solver.radial_eigenvalues, **fields}
    path = directory/(case.name+'.npz')
    np.savez_compressed(path, **arrays)
    modal_energy = solver.energies(state)
    row = {'name': case.name, 'completed': True, 'configuration': case.configuration(), 'provenance': source,
           'elapsed_seconds': perf_counter()-start, 'end_time_s': end, 'dt_s': dt,
           'continuum_frequency_real_s-1': float(reference.omega.real), 'continuum_damping_s-1': float(-reference.omega.imag),
           'phase_frequency_relative_error': float(phase_error), 'measured_damping_s-1': float(damping),
           'damping_error_over_frequency': float((damping+reference.omega.imag)/reference.omega.real),
           'initial_errors': initial_errors, 'final_errors': errors, 'time_relative_energy_norm_error': temporal_error,
           'reference_heat_j': exact_heat, 'heat_error_over_initial_energy': (heat-exact_heat)/initial_energy,
           'physical_energy_j': physical, 'physical_vs_modal_relative_error': max(abs(physical[k]-modal_energy[k]) for k in physical)/initial_energy,
           'physical_quadrature_relative_change': max(abs(physical[k]-physical12[k]) for k in physical)/initial_energy,
           'error_quadrature_absolute_change': max(abs(errors[k]-errors12[k]) for k in errors),
           'heat_quadrature_change_over_energy': abs(exact_heat-exact_heat12)/initial_energy,
           'history': history, 'fields': {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                                         'units': 'SI; magnetic T, velocity m/s, current A/m², volume_weights m³, modes documented in hall_modes'},
           'boundary': 'zero a, Btheta, velocity stream and vtheta at r=0,R; periodic z; no vacuum coupling'}
    if source != provenance():
        raise RuntimeError('Hall source changed during measurement')
    write_record(directory/(case.name+'.json'), row)
    print('Completed', case.name, flush=True)
    return row


def run(directory='outputs/hall_verification'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    rows = [measure(case, directory) for case in protocol()]
    if provenance() != source:
        raise RuntimeError('Source changed across Hall cases')
    write_record(directory/'measurements.json', {'completed': True, 'provenance': source, 'cases': rows})


if __name__ == '__main__':
    run()
