"""Actual isolated pressure/enthalpy energy evolution on primal and dual states."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_mass import DualMass
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_midpoint import step
from plasma_bench.electron_energy import energy_rate
from plasma_bench.thermal_means import logarithmic_mean
from .electron_reference import ElectronReference


def measure(case):
    ref = ElectronReference(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    mapping = DualMass(mesh); geo = mapping.geometry
    f = ref.fields(geo.r, geo.z)
    model = ElectronPressure(mesh, f['density'], gamma=ref.gamma)
    u, g = model.from_pressure(f['pressure']), ref.cell_average(mesh, 'g')
    pi0 = case['ion_pressure']
    total = pi0/(ref.gamma-1)+mapping.primal_average(model.energy_weight*u)+.5*geo.magnetic_mass*g*g/mesh.volumes
    arrays = dict(initial_g=g.copy(), initial_entropy=u.copy(), initial_total_energy=total.copy())
    def diagnostics(time, iterations=0, residual=0.):
        magnetic = .5*geo.magnetic_mass*g*g/mesh.volumes
        electron = mapping.primal_average(model.energy_weight*u)
        ion = total-magnetic-electron
        if np.any(ion <= 0):
            raise RuntimeError('Isolated pressure step produced nonpositive ion energy; no repair')
        return dict(time_s=time, magnetic_energy_j=float(np.sum(mesh.volumes*magnetic)), electron_internal_energy_j=float(np.sum(mesh.volumes*electron)),
                    ion_internal_energy_j=float(np.sum(mesh.volumes*ion)), total_energy_j=float(np.sum(mesh.volumes*total)),
                    entropy_proxy_j=float(np.sum(geo.volumes*u)), entropy_log_particles=float(np.sum(geo.volumes*model.density*np.log(u/(model.density*E)))),
                    min_ion_pressure_pa=float((ref.gamma-1)*ion.min()), min_electron_pressure_pa=float(model.pressure(u).min()),
                    nonlinear_iterations=iterations, nonlinear_residual=residual)
    initial = diagnostics(0.); history = {k: [v] for k, v in initial.items()}
    dt = case['end_time']/case['steps']
    for i in range(case['steps']):
        result = step(model, g, u, dt, tolerance=case['tolerance'], max_iterations=case['pressure_iterations'])
        total += dt*energy_rate(model, .5*(g+result.g), logarithmic_mean(u, result.entropy_density))
        g, u = result.g, result.entropy_density
        last = diagnostics((i+1)*dt, result.iterations, result.residual)
        for k, v in last.items():
            history[k].append(v)
    pi = (ref.gamma-1)*(total-mapping.primal_average(model.energy_weight*u)-.5*geo.magnetic_mass*g*g/mesh.volumes)
    metrics = dict(ion_pressure_rms_change_pa=float(np.sqrt(np.sum(mesh.volumes*(pi-pi0)**2)/np.sum(mesh.volumes))),
                   ion_pressure_max_change_pa=float(np.max(abs(pi-pi0))),
                   total_energy_relative_change=(last['total_energy_j']-initial['total_energy_j'])/initial['total_energy_j'])
    arrays.update(final_g=g, final_entropy=u, final_total_energy=total, final_ion_pressure=pi, number_density=model.density,
                  primal_volumes=mesh.volumes, dual_volumes=geo.volumes, radial_faces=mesh.rf, axial_faces=mesh.zf)
    return dict(initial=initial, final=last, metrics=metrics, accepted_steps=case['steps'], end_time_s=case['end_time'], history_columns=history), arrays
