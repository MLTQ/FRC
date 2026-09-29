"""Actual coupled warm Hall trajectories and complete accepted histories."""
import numpy as np
from .warm_hall_initial import initialize
from .warm_hall_norms import wave_errors


def measure(case):
    reference, mesh, model = initialize(case)
    initial_state, initial_u = model.state.copy(), model.electron_entropy.copy()
    initial_pe = model.thermal.pressure(initial_state, initial_u)
    wave = case['kind'] == 'wave'
    initial_errors = wave_errors(model, reference)[0] if wave else None
    end = case['phase']/abs(reference.omega) if wave else case['end_time']
    cap = case['phase_step']/abs(reference.omega) if wave else case['dt_max']
    keys = ('time_s', 'steps', 'mass_kg', 'total_energy_j', 'internal_energy_j', 'kinetic_energy_j', 'magnetic_energy_j',
            'electron_internal_energy_j', 'ion_internal_energy_j', 'min_density_kg_m3', 'min_pressure_pa',
            'min_electron_pressure_pa', 'min_ion_pressure_pa', 'max_face_divergence_t_m', 'physical_ohmic_heat_j',
            'electron_entropy_proxy_integral_j', 'electron_log_entropy_particles', 'hall_nonlinear_iterations',
            'pressure_nonlinear_iterations', 'max_hall_nonlinear_residual', 'max_hall_linear_residual', 'max_pressure_nonlinear_residual')
    initial = model.diagnostics(); history = {k: [initial[k]] for k in keys}
    def record():
        d = model.diagnostics()
        for key in keys:
            history[key].append(d[key])
    if case['steps'] is not None:
        for _ in range(case['steps']):
            model.step(end/case['steps']); record()
    else:
        while model.time < end:
            remaining = end-model.time
            if remaining <= 4*np.finfo(float).eps*end:
                break
            model.step(min(remaining, cap, model.timestep_limit())); record()
            if model.steps % 250 == 0:
                print('Progress', case['name'], model.steps, model.time/end, flush=True)
    pe, pi = model.thermal.validate(model.state, model.electron_entropy)
    errors = wave_errors(model, reference)[0] if wave else None
    arrays = dict(initial_conserved=initial_state.conserved, final_conserved=model.state.conserved,
                  initial_psi=initial_state.psi, final_psi=model.state.psi, initial_entropy=initial_u,
                  final_entropy=model.electron_entropy, initial_electron_pressure=initial_pe, final_electron_pressure=pe,
                  final_ion_pressure=pi, primal_volumes=mesh.volumes, dual_volumes=model.thermal.mapping.volumes,
                  radial_faces=mesh.rf, axial_faces=mesh.zf)
    return dict(initial=initial, final=model.diagnostics(), initial_errors=initial_errors, metrics=errors,
                end_time_s=end, actual_time_s=model.time, accepted_steps=model.steps, history_columns=history), arrays
