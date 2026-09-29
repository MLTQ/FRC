"""Coupled Hall-fluid evolution, physical wave errors and conservative histories."""
from time import perf_counter
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.mhd import MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.hall_fluid import HallMHD
from plasma_bench.toroidal_state import primitive
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction
from .compressible_hall_reference import CompressibleHallWave
from .hall_manufactured import HallManufactured
from .compatible_hall_references import wave_state, nonlinear_state
from .compatible_hall_records import provenance, relative, save


def wave_errors(solver, reference, time):
    geo, mesh = solver.geometry, solver.mesh
    q = primitive(mesh, solver.state, solver.parameters.gamma)
    rec = ToroidalReconstruction(mesh, q, order=solver.parameters.spatial_order)
    actual = np.empty(geo.r.shape+(8,))
    for i in range(geo.order):
        radius = geo.r[:, 0, i, 0]
        for j in range(geo.order):
            offset = (geo.z[0, 0, 0, j]-mesh.z[0])/mesh.dz
            actual[:, :, i, j] = rec.primitive(rec.raw(radius, offset), rec.auxiliary(radius, offset), radius[:, None])
    expected = reference.fields(geo.r, geo.z, time)
    velocity = np.stack((actual[..., 1], geo.r*actual[..., 6], actual[..., 2]), axis=-1)
    magnetic = geo.magnetic_samples(solver.state.psi, solver.state.conserved[..., 5])-np.array([0., 0., reference.bias])
    result = dict(magnetic_relative_L2=relative(magnetic-expected['magnetic'], expected['magnetic'], geo.weight[..., None]),
                  velocity_relative_L2=relative(velocity-expected['velocity'], expected['velocity'], geo.weight[..., None]),
                  density_relative_perturbation_L2=relative(actual[..., 0]-reference.density-expected['density'], expected['density'], geo.weight),
                  pressure_relative_perturbation_L2=relative(actual[..., 3]-reference.pressure-expected['pressure'], expected['pressure'], geo.weight))
    exact = wave_state(mesh, reference, time)
    result['g_cell_relative_L2'] = relative(solver.state.conserved[..., 5]-exact.conserved[..., 5], exact.conserved[..., 5], mesh.volumes)
    return result, exact


def measure(case, directory):
    source, start = provenance(), perf_counter()
    is_wave = case['kind'] == 'wave'
    ref = (CompressibleHallWave if is_wave else HallManufactured)(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    mi = ref.ion_mass if is_wave else 39.948*atomic_mass
    gamma = ref.gamma if is_wave else 5/3
    initial = wave_state(mesh, ref) if is_wave else nonlinear_state(mesh, ref, mi)
    solver = HallMHD(mesh, initial, mi, MHDParameters(gamma=gamma, resistivity=case['resistivity']),
                    hall=ref.hall if is_wave else True, quadrature_order=case['quadrature'], nonlinear_tolerance=case['nonlinear_tolerance'])
    end = case['phase']/abs(ref.omega) if is_wave else case['end_time']
    initial_errors = wave_errors(solver, ref, 0.)[0] if is_wave else None
    keys = ('time_s', 'mass_kg', 'total_energy_j', 'internal_energy_j', 'kinetic_energy_j', 'magnetic_energy_j',
            'angular_momentum_kg_m2_s', 'physical_ohmic_heat_j', 'min_density_kg_m3', 'min_pressure_pa',
            'max_face_divergence_t_m', 'hall_nonlinear_iterations', 'max_hall_nonlinear_residual', 'max_hall_linear_residual')
    history = {key: [] for key in keys}
    def record():
        out = solver.diagnostics()
        for key in keys:
            history[key].append(out[key])
    record()
    while solver.time < end:
        remaining = end-solver.time
        if remaining <= 4*np.finfo(float).eps*end:
            solver.time = end
            break
        solver.step(min(remaining, solver.timestep_limit(), case['dt_max'] if case['dt_max'] is not None else np.inf))
        record()
        if solver.steps % 500 == 0:
            print('Progress', case['name'], solver.steps, 'steps', solver.time/end, flush=True)
    errors, exact = wave_errors(solver, ref, end) if is_wave else (None, None)
    arrays = dict(initial_conserved=initial.conserved, initial_psi=initial.psi,
                  final_conserved=solver.state.conserved, final_psi=solver.state.psi,
                  radial_faces=mesh.rf, axial_faces=mesh.zf, volumes=mesh.volumes)
    if exact is not None:
        arrays.update(reference_conserved=exact.conserved, reference_psi=exact.psi)
    result = save(directory, case, dict(kind='coupled_evolution', end_time_s=end, elapsed_seconds=perf_counter()-start,
                   initial=solver.initial_integrals, final=solver.diagnostics(), initial_errors=initial_errors,
                   errors=errors, history_columns=history), arrays, source)
    print('Completed', case['name'], solver.steps, 'steps', errors, flush=True)
    return result
