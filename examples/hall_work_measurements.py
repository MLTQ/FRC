"""Reproducible spatial Hall current, induction and local/global work audits."""
from pathlib import Path
import hashlib
import json
import sys
from time import perf_counter
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.hall_curl import HallCurl
from plasma_bench.hall_induction import hall_induction
from plasma_bench.hall_poynting import local_energy_audit
from plasma_bench.toroidal_state import radius_squared
from .compressible_hall_reference import CompressibleHallWave
from .hall_manufactured import HallManufactured
from .hall_work_snapshots import snapshot, wave_snapshot, wave_hall_rate
from .hall_work_cases import protocol
from .verification_records import write_record


def provenance():
    names = [f'plasma_bench/{name}.py' for name in ('hall_curl', 'hall_induction', 'hall_poynting', 'mhd_mesh', 'toroidal_state', 'mhd_state', 'magnetic_quadrature')]
    names += [f'examples/{name}.py' for name in ('compressible_hall_reference', 'hall_manufactured', 'hall_work_snapshots', 'hall_work_cases', 'hall_work_measurements', 'verification_records')]
    files = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}
    return dict(source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                source_files=files, python=sys.version, numpy=np.__version__, scipy=scipy.__version__)


def relative(a, b, weight):
    base = np.sum(weight*b*b)
    return float(np.sqrt(np.sum(weight*(a-b)**2)/base)) if base > 0 else None


def measure(case, directory):
    start, source = perf_counter(), provenance()
    kind = case['kind']
    reference = (HallManufactured if kind == 'nonlinear' else CompressibleHallWave)(**case['reference'])
    mesh = MHDMesh(reference.radius, reference.length/2, case['nr'], case['nz'])
    geometry = HallCurl(mesh, case['quadrature'])
    make = snapshot if kind == 'nonlinear' else wave_snapshot
    initial, control = make(reference, mesh, 10), make(reference, mesh, 14)
    exact = reference.fields(geometry.r, geometry.z)
    if kind == 'wave':
        exact['induction'], _ = wave_hall_rate(reference, geometry.r, geometry.z)
        exact['number_density'] = np.full(geometry.r.shape, reference.number_density)
    density = exact['number_density'] if case['density_sampling'] == 'quadrature' else initial['density']
    result = hall_induction(geometry, initial['psi'], initial['g'], density)
    actual_rate = geometry.magnetic_samples(result.psi_rate, result.g_rate)
    magnetic, total, sr, sz = local_energy_audit(geometry, initial['psi'], initial['g'], result)
    weight = geometry.weight[..., None]
    magnetic_scale = initial['magnetic_dot_induction']/MU0
    energy_scale2 = np.sum(mesh.volumes*magnetic_scale**2)
    energy_defect2 = np.sum(mesh.volumes*(magnetic-total)**2)
    projected = geometry.sample_electric(result.electric_coefficients)
    raw_work = np.sum(result.current_samples*result.electric_samples, axis=-1)
    projected_work = np.sum(result.current_samples*projected, axis=-1)
    metrics = dict(current_relative_L2=relative(result.current_samples, exact['current'], weight),
                   induction_relative_L2=relative(actual_rate, exact['induction'], weight),
                   g_rate_cell_relative_L2=relative(result.g_rate, initial['g_rate'], mesh.volumes),
                   toroidal_rate_cell_relative_L2=relative(result.g_rate, initial['g_rate'], mesh.volumes*radius_squared(mesh)),
                   axis_g_rate_relative_L2=relative(result.g_rate[0], initial['g_rate'][0], np.ones(mesh.nz)),
                   psi_rate_max_relative=float(np.max(abs(result.psi_rate-initial['psi_rate']))/np.max(abs(initial['psi_rate']))) if np.max(abs(initial['psi_rate'])) > 0 else None,
                   local_energy_defect_relative_L2=float(np.sqrt(energy_defect2/energy_scale2)) if energy_scale2 > 0 else None,
                   local_energy_defect_rms_w_m3=float(np.sqrt(energy_defect2/np.sum(mesh.volumes))),
                   global_magnetic_work_w=result.magnetic_work_w, absolute_work_scale_w=result.absolute_work_scale_w,
                   global_magnetic_work_relative=abs(result.magnetic_work_w)/result.absolute_work_scale_w,
                   total_energy_flux_integral_w=float(np.sum(mesh.volumes*total)),
                   raw_pointwise_work_max_w_m3=float(np.max(abs(raw_work))),
                   projected_pointwise_work_rms_w_m3=float(np.sqrt(np.sum(geometry.weight*projected_work**2)/np.sum(geometry.weight))),
                   poloidal_energy_rate_w=float(np.sum(mesh.volumes*(magnetic-radius_squared(mesh)*initial['g']*result.g_rate/MU0))),
                   toroidal_energy_rate_w=float(np.sum(mesh.volumes*radius_squared(mesh)*initial['g']*result.g_rate/MU0)),
                   axis_density_min_m3=float(exact['number_density'][0].min()),
                   cell_reference_quadrature_relative_change=max(relative(initial[key], control[key], np.ones_like(initial[key])) or 0. for key in ('g', 'g_rate', 'density')))
    if kind == 'nonlinear':
        metrics['electric_raw_relative_L2'] = relative(result.electric_samples, exact['electric'], weight)
        metrics['electric_projected_relative_L2'] = relative(projected, exact['electric'], weight)
    # Both signs separate the nonlinear quadratic Hall contribution from its
    # linearization about a uniform bias without subtracting a fitted reference.
    rc, _ = mesh.corners
    bias = .5*reference.bias*rc*rc
    negative = hall_induction(geometry, 2*bias-initial['psi'], -initial['g'], density)
    nonlinear = geometry.magnetic_samples(.5*(result.psi_rate+negative.psi_rate), .5*(result.g_rate+negative.g_rate))
    # Independently exploit bilinearity: Hall(b,b) evaluates the quadratic
    # contribution directly, avoiding cancellation of O(bias*perturbation).
    pure = hall_induction(geometry, initial['psi']-bias, initial['g'], density)
    direct = geometry.magnetic_samples(pure.psi_rate, pure.g_rate)
    metrics['quadratic_induction_rms_t_s'] = float(np.sqrt(np.sum(weight*direct**2)/np.sum(geometry.weight)))
    metrics['quadratic_symmetric_difference_rms_t_s'] = float(np.sqrt(np.sum(weight*nonlinear**2)/np.sum(geometry.weight)))
    metrics['quadratic_symmetric_difference_relative_error'] = relative(nonlinear, direct, weight)
    arrays = dict(psi=initial['psi'], g=initial['g'], psi_rate=result.psi_rate, g_rate=result.g_rate,
                  reference_psi_rate=initial['psi_rate'], reference_g_rate=initial['g_rate'],
                  density_quadrature=np.broadcast_to(density[..., None, None] if density.shape == initial['g'].shape else density, geometry.r.shape),
                  magnetic_energy_rate=magnetic, total_energy_flux_rate=total, radial_poynting=sr, axial_poynting=sz,
                  electric_coefficients=result.electric_coefficients, radial_faces=mesh.rf, axial_faces=mesh.zf)
    path = directory/(case['name']+'.npz')
    np.savez_compressed(path, **arrays)
    record = dict(completed=True, kind='instantaneous_spatial_audit', configuration=case, provenance=source,
                  metrics=metrics, elapsed_seconds=perf_counter()-start,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; psi Wb/(2pi), g T/m, rates per second; no evolved state'))
    if provenance() != source:
        raise RuntimeError('Source changed during Hall spatial measurement')
    write_record(directory/(case['name']+'.json'), record)
    print('Completed', case['name'], 'induction error', metrics['induction_relative_L2'], flush=True)
    return record


def run(directory='outputs/hall_work_verification'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    records = [measure(case, directory) for case in protocol()]
    if provenance() != source:
        raise RuntimeError('Source changed across Hall spatial study')
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, cases=records))


if __name__ == '__main__':
    run()
