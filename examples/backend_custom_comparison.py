"""Existing linear Hall component on the external-engine continuum benchmark."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import proton_mass
from plasma_bench.hall_parameters import HallParameters
from plasma_bench.hall_modes import LinearHallModes
from .backend_wave_reference import ParallelWave, relative_error


def run(directory='outputs/engine_comparison'):
    rows = []
    for h in (-1, 1):
        wave = ParallelWave(helicity=h)
        parameters = HallParameters(wave.density, proton_mass, wave.bias, 0., 'hall')
        solver = LinearHallModes([wave.k], [0.], parameters)
        v_ratio = -wave.k*parameters.alfven/wave.omega
        initial = wave.amplitude*np.array([[[1., h, v_ratio, h*v_ratio]]], complex)
        stop = 2/wave.omega
        for dt_omega in (.005, .0025):
            count = int(np.ceil(stop*wave.omega_i/dt_omega))
            state = initial.copy()
            e0 = sum(solver.energies(state).values())
            begin = time.perf_counter()
            for _ in range(count):
                state = solver.step(state, stop/count).state
            elapsed = time.perf_counter()-begin
            exact = initial*np.exp(-1j*wave.omega*stop)
            rows.append({'helicity': h, 'steps': count, 'stop_s': stop,
                'dt_s': stop/count, 'dt_omega': dt_omega,
                'relative_state_error': relative_error(state, exact),
                'relative_energy_change': (sum(solver.energies(state).values())-e0)/e0,
                'evolution_s': elapsed,
                'final_real': state.real.tolist(), 'final_imag': state.imag.tolist()})
    names = ['plasma_bench/hall_modes.py', 'plasma_bench/hall_parameters.py',
             'examples/backend_wave_reference.py', 'examples/backend_custom_comparison.py']
    result = {'completed': True, 'scope': 'Existing linear four-mode Hall component; no spatial grid or particles. Runtime is not comparable to a PIC evolution.',
              'cases': rows, 'sources': {n: hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}}
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    (path/'custom_hall.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    run()
