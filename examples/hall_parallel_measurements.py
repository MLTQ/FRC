"""Independent parallel circular-wave dispersion/damping and time convergence."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import proton_mass as MP, elementary_charge as E, mu_0 as MU0
from plasma_bench.hall_parameters import HallParameters
from plasma_bench.hall_modes import LinearHallModes
from .hall_wave_reference import frequency
from .verification_records import write_record


def run(directory='outputs/hall_verification'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    names = ['plasma_bench/hall_parameters.py', 'plasma_bench/hall_modes.py',
             'examples/hall_wave_reference.py', 'examples/hall_parallel_measurements.py', 'examples/verification_records.py']
    hashes = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}
    rows = []
    for model in ('hall', 'two_fluid'):
        p = HallParameters(1e18, MP, .001, 3e-5, model)
        di = np.sqrt(MP/(MU0*p.number_density*E**2))
        for kd in (.1, 1., 10., 1000.):
            k = kd/di
            solver = LinearHallModes([k], [0.], p)
            for h in (-1, 1):
                omega = frequency(p.number_density, MP, p.bias, p.resistivity, k, 0., h, model)
                ratio = -k*p.alfven/omega
                initial = np.array([[[1., h, ratio, h*ratio]]], complex)*1e-8
                time = 2/omega.real
                e0 = sum(solver.energies(initial, .5).values())
                errors, refinements = [], []
                for count in (40, 80, 160):
                    state, heat = initial.copy(), 0.
                    for _ in range(count):
                        step = solver.step(state, time/count, .5)
                        state, heat = step.state, heat+step.heat
                    exact = initial*np.exp(-1j*omega*time)
                    amplitude = state[0, 0, 0]/initial[0, 0, 0]
                    error = float(np.linalg.norm(state-exact)/np.linalg.norm(exact))
                    errors.append(error)
                    refinements.append({'steps': count, 'dt_s': time/count, 'relative_state_error': error,
                                        'phase_frequency_relative_error': float(-np.angle(amplitude*np.exp(1j*omega.real*time))/2),
                                        'damping_error_over_frequency': float((-np.log(abs(amplitude))/time+omega.imag)/omega.real),
                                        'relative_energy_balance': float((sum(solver.energies(state, .5).values())+heat-e0)/e0),
                                        'heat_error_over_initial_energy': float(heat/e0+np.expm1(2*omega.imag*time)),
                                        'final_coefficients_real': state.real.tolist(), 'final_coefficients_imag': state.imag.tolist()})
                rows.append({'model': model, 'helicity': h, 'k_di': kd, 'k_m-1': float(k),
                             'number_density_m-3': p.number_density, 'ion_mass_kg': p.ion_mass,
                             'bias_t': p.bias, 'resistivity_ohm_m': p.resistivity,
                             'frequency_real_s-1': float(omega.real), 'damping_s-1': float(-omega.imag),
                             'time_error_ratios': [errors[0]/errors[1], errors[1]/errors[2]], 'refinements': refinements})
    if any(hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha for name, sha in hashes.items()):
        raise RuntimeError('Source changed during parallel-wave measurement')
    result = {'completed': True, 'scope': 'Cartesian parallel circular transverse waves; not cylindrical formation',
              'energy_units': 'J/m³ averaged over the real sinusoid, with drift energy',
              'source_files': hashes, 'source_sha256': hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(), 'cases': rows}
    write_record(directory/'parallel.json', result)
    print('Completed', len(rows), 'parallel modes, three time steps each')
    return result


if __name__ == '__main__':
    run()
