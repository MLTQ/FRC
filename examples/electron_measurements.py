"""Run independent electron-pressure spatial and actual evolution measurements."""
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_current import DualCurrent
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_midpoint import step, ElectronConvergenceError
from .electron_reference import ElectronReference, drift_pressure
from .electron_cases import protocol
from .electron_records import provenance, relative, save
from .verification_records import write_record


def initialize(case):
    ref = ElectronReference(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    geo = DualCurrent(mesh)
    fields = ref.fields(geo.r, geo.z)
    if case['kind'] == 'drift':
        n = np.full_like(geo.r, ref.density)
        g = np.full((mesh.nr, mesh.nz), case['slope'])
        p = drift_pressure(geo, 0., amplitude=case['amplitude'], slope=case['slope'])
        wall = case['slope']*mesh.radius
    else:
        n, p = fields['density'], fields['pressure']
        g = ref.cell_average(mesh, 'g')
        if case['kind'] == 'battery':
            g = np.zeros_like(g)
        wall = 0.
    model = ElectronPressure(mesh, n, gamma=ref.gamma, density_scale=case['density_scale'], wall_btheta=wall)
    return ref, mesh, model, g, model.from_pressure(p)


def spatial(case, ref, mesh, model, g, u):
    geo = model.geometry
    dg, du, er, ez = model.rhs(g, u)
    expected = ref.cell_average(mesh, 'g_rate')
    check = ref.cell_average(mesh, 'g_rate', 14)
    nodal = ref.fields(geo.r, geo.z)
    ir, iz = geo.current_flux(g, 0.)
    x = geo.dual_x_faces
    # Independent continuum line/annular current integrals use analytic r Btheta.
    rr = np.sqrt(2*x)[:, None]
    true_ir = ref.fields(rr, (mesh.zf[:-1]-.5*mesh.dz)[None])['stream']-ref.fields(rr, (mesh.zf[:-1]+.5*mesh.dz)[None])['stream']
    true_iz = ref.fields(np.sqrt(2*x[1:])[:, None], mesh.z[None])['stream']-ref.fields(np.sqrt(2*x[:-1])[:, None], mesh.z[None])['stream']
    summation = max(float(abs(ir).max()), float(abs(iz).max()), np.finfo(float).tiny)
    divergence = np.diff(ir, axis=0)+iz-np.roll(iz, 1, axis=1)
    thermal_rate = geo.volumes*model.energy_weight*du
    magnetic_rate = geo.magnetic_mass*g*dg
    scale = float(np.sum(abs(thermal_rate))+np.sum(abs(magnetic_rate)))
    metrics = dict(battery_relative_L2=relative(dg-expected, expected, mesh.volumes),
                   battery_rms_error_t_m_s=float(np.sqrt(np.sum(mesh.volumes*(dg-expected)**2)/np.sum(mesh.volumes))),
                   axis_battery_relative_L2=relative(dg[0]-expected[0], expected[0], np.ones(mesh.nz)),
                   thermal_rate_relative_L2=relative(du-nodal['entropy_rate']* (ref.density/case['density_scale'])**(1-ref.gamma),
                                                     nodal['entropy_rate']*(ref.density/case['density_scale'])**(1-ref.gamma), geo.volumes),
                   radial_current_flux_relative_L2=relative(ir-true_ir, true_ir, np.ones_like(ir)),
                   axial_current_flux_relative_L2=relative(iz-true_iz, true_iz, np.ones_like(iz)),
                   axis_axial_current_relative_L2=relative(iz[0]-true_iz[0], true_iz[0], np.ones(mesh.nz)),
                   charge_cancellation_relative=float(abs(divergence).max()/summation),
                   work_cancellation_relative=float(abs(np.sum(thermal_rate)+np.sum(magnetic_rate))/scale) if scale else 0.,
                   logarithmic_entropy_rate=float(np.sum(geo.volumes*model.density/u*du)),
                   reference_quadrature_change=relative(check-expected, check, mesh.volumes))
    return dict(kind='spatial', metrics=metrics), dict(g=g, entropy_density=u, battery_rate=dg, thermal_rate=du,
                   reference_battery_rate=expected, reference_thermal_rate=nodal['entropy_rate'], current_radial=ir, current_axial=iz,
                   reference_current_radial=true_ir, reference_current_axial=true_iz, electric_radial_coefficient=er, electric_axial=ez)


def evolution(case, ref, mesh, model, g, u):
    geo = model.geometry
    initial_g, initial_u = g.copy(), u.copy()
    initial = model.diagnostics(g, u)
    end = ref.length/(2*abs(case['slope'])/(MU0*ref.density*E)) if case['kind'] == 'drift' else case['end_time']
    dt, time, iterations, residual = end/case['steps'], 0., 0, 0.
    history = {key: [value] for key, value in initial.items()}
    history['time_s'] = [0.]
    failure = None
    for i in range(case['steps']):
        before_g, before_u = g.copy(), u.copy()
        try:
            result = step(model, g, u, dt, tolerance=case['tolerance'])
        except ElectronConvergenceError as error:
            np.testing.assert_array_equal(g, before_g); np.testing.assert_array_equal(u, before_u)
            failure = str(error)
            break
        g, u = result.g, result.entropy_density
        iterations += result.iterations; residual = max(residual, result.residual)
        time = end*(i+1)/case['steps']
        for key, value in model.diagnostics(g, u).items():
            history[key].append(value)
        history['time_s'].append(time)
    if (failure is not None) != case['expected_failure']:
        raise RuntimeError(f'Unexpected completion/failure in {case["name"]}: {failure}')
    metrics = {}
    if case['kind'] == 'drift':
        exact = drift_pressure(geo, time, amplitude=case['amplitude'], slope=case['slope'])
        metrics['pressure_perturbation_relative_L2'] = relative(model.pressure(u)-exact, exact-80., geo.volumes)
        metrics['magnetic_relative_change'] = relative(g-initial_g, initial_g, mesh.volumes)
    particles = float(np.sum(geo.volumes*model.density))
    final = model.diagnostics(g, u)
    metrics['logarithmic_entropy_change_per_particle'] = (final['entropy_log_particle_units']-initial['entropy_log_particle_units'])/particles
    metrics['total_energy_relative_change'] = (final['total_energy_j']-initial['total_energy_j'])/initial['total_energy_j']
    return dict(kind='evolution', initial=initial, final=final, end_time_s=end, actual_time_s=time,
                accepted_steps=len(history['time_s'])-1, iterations=iterations, max_nonlinear_residual=residual,
                failure=failure, metrics=metrics, history_columns=history), dict(initial_g=initial_g, final_g=g,
                initial_entropy_density=initial_u, final_entropy_density=u, initial_pressure=model.pressure(initial_u), final_pressure=model.pressure(u))


def run(directory='outputs/electron_pressure_verification'):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source, records = provenance(), []
    for case in protocol():
        start = perf_counter()
        ref, mesh, model, g, u = initialize(case)
        result, arrays = (spatial if case['kind'] == 'spatial' else evolution)(case, ref, mesh, model, g, u)
        arrays.update(number_density=model.density, dual_volumes=model.geometry.volumes, primal_volumes=mesh.volumes,
                      radial_faces=mesh.rf, axial_faces=mesh.zf)
        records.append(save(directory, case, source, result | dict(elapsed_seconds=perf_counter()-start), arrays))
        print(case['name'], result.get('failure') or result['metrics'], flush=True)
    if provenance() != source:
        raise RuntimeError('Sources changed during electron protocol')
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, cases=records))


if __name__ == '__main__':
    run()
