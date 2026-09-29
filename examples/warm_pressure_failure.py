"""Preserve accepted and rejected pressure-only states at the newly failing mesh."""
from copy import deepcopy
from pathlib import Path
import hashlib
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_mass import DualMass
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_midpoint import step
from plasma_bench.electron_energy import energy_rate
from plasma_bench.thermal_means import logarithmic_mean
from .electron_reference import ElectronReference
from .warm_hall_controls import protocol, source_hash as control_hash
from .warm_hall_records import provenance
from .verification_records import write_record


def run(count, directory='outputs/warm_pressure_failure'):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = provenance(); driver = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    original = next(c for c in protocol() if c['name'] == 'pressure_fine_mesh_128')
    case = deepcopy(original); case.update(name=f'pressure_mesh128_time{count}', steps=count)
    ref = ElectronReference(**case['reference']); mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    mapping = DualMass(mesh); geo = mapping.geometry; f = ref.fields(geo.r, geo.z)
    model = ElectronPressure(mesh, f['density'], gamma=ref.gamma)
    u, g = model.from_pressure(f['pressure']), ref.cell_average(mesh, 'g')
    total = case['ion_pressure']/(ref.gamma-1)+mapping.primal_average(model.energy_weight*u)+.5*geo.magnetic_mass*g*g/mesh.volumes
    def diagnostics(g, u, total, time, iterations=0, residual=0.):
        magnetic = .5*geo.magnetic_mass*g*g/mesh.volumes
        electron = mapping.primal_average(model.energy_weight*u)
        pi = (ref.gamma-1)*(total-magnetic-electron)
        idx = np.unravel_index(np.argmin(pi), pi.shape)
        return dict(time_s=time, magnetic_energy_j=float(np.sum(mesh.volumes*magnetic)), electron_internal_energy_j=float(np.sum(mesh.volumes*electron)),
                    ion_internal_energy_j=float(np.sum(mesh.volumes*pi)/(ref.gamma-1)), total_energy_j=float(np.sum(mesh.volumes*total)),
                    entropy_proxy_j=float(np.sum(geo.volumes*u)), entropy_log_particles=float(np.sum(geo.volumes*model.density*np.log(u/(model.density*E)))),
                    min_ion_pressure_pa=float(pi[idx]), min_electron_pressure_pa=float(model.pressure(u).min()),
                    minimum_ion_cell=[int(i) for i in idx], nonlinear_iterations=iterations, nonlinear_residual=residual), pi
    initial, _ = diagnostics(g, u, total, 0.)
    history = {k: [v] for k, v in initial.items() if k != 'minimum_ion_cell'}
    arrays = dict(initial_g=g.copy(), initial_entropy=u.copy(), initial_total_energy=total.copy())
    dt = case['end_time']/count; rejected = None; accepted = 0; last = initial
    for i in range(count):
        result = step(model, g, u, dt, tolerance=case['tolerance'], max_iterations=case['pressure_iterations'])
        trial_total = total+dt*energy_rate(model, .5*(g+result.g), logarithmic_mean(u, result.entropy_density))
        trial, pi = diagnostics(result.g, result.entropy_density, trial_total, (i+1)*dt, result.iterations, result.residual)
        if trial['min_ion_pressure_pa'] <= 0:
            rejected = trial | dict(reason='nonpositive ion energy from conservative local flux; no floor or reset')
            arrays.update(rejected_g=result.g, rejected_entropy=result.entropy_density, rejected_total_energy=trial_total, rejected_ion_pressure=pi)
            break
        g, u, total, last, accepted = result.g, result.entropy_density, trial_total, trial, i+1
        for k in history:
            history[k].append(last[k])
    last, pi = diagnostics(g, u, total, accepted*dt, last['nonlinear_iterations'], last['nonlinear_residual'])
    arrays.update(final_g=g, final_entropy=u, final_total_energy=total, final_ion_pressure=pi, number_density=model.density,
                  primal_volumes=mesh.volumes, dual_volumes=geo.volumes, radial_faces=mesh.rf, axial_faces=mesh.zf)
    assert provenance() == source and hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == driver
    path = directory/(case['name']+'.npz'); np.savez_compressed(path, **arrays)
    record = dict(completed=accepted == count, experiment_complete=True, configuration=case, original_control_configuration=original,
                  provenance=source, driver_sha256=driver, control_source_sha256=control_hash(), initial=initial, final=last,
                  accepted_steps=accepted, end_time_s=case['end_time'], actual_time_s=accepted*dt, history_columns=history, rejected_trial=rejected,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; accepted and rejected arrays separately labeled'))
    write_record(directory/(case['name']+'.json'), record)
    print(case['name'], 'accepted', accepted, 'last time', accepted*dt, 'rejected', rejected, flush=True)


if __name__ == '__main__':
    for count in (320, 640, 1280):
        run(count)
