"""Independent Bessel decay and slow-wave amplitude controls for Hall integration."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import numpy as np
from scipy.constants import atomic_mass, mu_0 as MU0
from scipy.special import j0, j1, jn_zeros
from plasma_bench.mhd import MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.hall_fluid import HallMHD
from plasma_bench.toroidal_state import from_primitive, radius_squared
from .compatible_hall_cases import protocol as main_protocol
from .compatible_hall_evolution import measure as evolve
from .compatible_hall_records import provenance, relative, save
from .verification_records import write_record


def protocol():
    source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    base = next(c for c in main_protocol() if c['name'] == 'wave_b0_n16')
    cases = [base | dict(name=f'slow_amplitude_{amplitude}', control_source_sha256=source,
                        reference=base['reference'] | dict(amplitude=amplitude)) for amplitude in (.1, .01)]
    for nr, steps in ((8, 64), (16, 64), (32, 64), (64, 64), (16, 32), (16, 128), (16, 256)):
        cases.append(dict(name=f'decay_n{nr}_s{steps}', kind='resistive_decay', nr=nr, nz=4, steps=steps,
                          radius=.035, length=.24, number_density=1e20, ion_mass=39.948*atomic_mass,
                          bias=.005, pressure=16., amplitude=1e-7, resistivity=1e-4, decay_exponent=.2,
                          control_source_sha256=source))
    return cases


def decay(case, directory):
    source, start = provenance(), perf_counter()
    mesh = MHDMesh(case['radius'], case['length']/2, case['nr'], case['nz'])
    q = jn_zeros(0, 1)[0]/mesh.radius
    rate = case['resistivity']/MU0*q*q
    duration = case['decay_exponent']/rate
    amplitude = case['amplitude']
    g0 = amplitude*(j0(q*mesh.rf[:-1])-j0(q*mesh.rf[1:]))[:, None]/(q*mesh.radial_measure[:, None])
    r, _ = mesh.corners
    initial = from_primitive(mesh, .5*case['bias']*r*r, case['number_density']*case['ion_mass'], case['pressure'],
                             toroidal_over_radius=g0, magnetic_quadrature='compatible')
    model = HallMHD(mesh, initial, case['ion_mass'], MHDParameters(resistivity=case['resistivity']), nonlinear_tolerance=1e-12)
    history = {key: [] for key in ('time_s', 'total_energy_j', 'physical_ohmic_heat_j', 'min_pressure_pa', 'max_hall_nonlinear_residual')}
    def record():
        value = model.diagnostics()
        for key in history:
            history[key].append(value[key])
    record()
    while model.time < duration:
        remaining = duration-model.time
        if remaining <= 4*np.finfo(float).eps*duration:
            model.time = duration
            break
        model.step(min(remaining, duration/case['steps'], model.timestep_limit()))
        record()
    exact = np.broadcast_to(g0*np.exp(-case['decay_exponent']), initial.conserved[..., 5].shape)
    continuum_energy = np.pi*case['length']*amplitude**2*mesh.radius**2*j1(q*mesh.radius)**2/(2*MU0)
    heat = continuum_energy*(-np.expm1(-2*case['decay_exponent']))
    out = model.diagnostics()
    metrics = dict(g_relative_L2=relative(model.state.conserved[..., 5]-exact, exact, mesh.volumes*radius_squared(mesh)),
                   axis_g_relative_error=float((model.state.conserved[0, 0, 5]-exact[0, 0])/exact[0, 0]),
                   physical_heat_relative_error=float((out['physical_ohmic_heat_j']-heat)/heat),
                   continuum_heat_j=float(heat), maximum_reference_field_t=float(amplitude),
                   actual_steps=model.steps)
    arrays = dict(initial_conserved=initial.conserved, initial_psi=initial.psi, final_conserved=model.state.conserved,
                  final_psi=model.state.psi, reference_g=exact, radial_faces=mesh.rf, axial_faces=mesh.zf, volumes=mesh.volumes)
    result = save(directory, case, dict(kind='resistive_decay', end_time_s=duration, elapsed_seconds=perf_counter()-start,
                  initial=model.initial_integrals, final=out, metrics=metrics, history_columns=history), arrays, source)
    print('Completed', case['name'], metrics, flush=True)
    return result


def run(directory='outputs/compatible_hall_controls'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    configurations = protocol()
    records = [(decay if c['kind'] == 'resistive_decay' else evolve)(c, directory) for c in configurations]
    if source != provenance() or configurations != protocol():
        raise RuntimeError('Engine or control protocol changed during Hall controls')
    write_record(directory/'measurements.json', dict(completed=True, provenance=source,
                 control_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), cases=records))


if __name__ == '__main__':
    run()
