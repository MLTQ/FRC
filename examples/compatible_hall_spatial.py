"""Independent compatible Hall current, axis and local transfer measurements."""
from time import perf_counter
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.hall_compatible import CompatibleHallCurl
from plasma_bench.hall_induction import hall_induction
from plasma_bench.compatible_poynting import energy_rate
from .hall_manufactured import HallManufactured
from .hall_work_snapshots import snapshot
from .compatible_hall_records import provenance, relative, save


def measure(case, directory):
    source, start = provenance(), perf_counter()
    ref = HallManufactured(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    geometry = CompatibleHallCurl(mesh, case['quadrature'])
    initial, control = snapshot(ref, mesh, 10), snapshot(ref, mesh, 14)
    exact = ref.fields(geometry.r, geometry.z)
    result = hall_induction(geometry, initial['psi'], initial['g'], exact['number_density'])
    rate = geometry.magnetic_samples(result.psi_rate, result.g_rate)
    magnetic = geometry.energy_rate(initial['psi'], initial['g'], result.psi_rate, result.g_rate)
    total = energy_rate(geometry, initial['psi'], initial['g'], result.electric_coefficients)
    local_reference = initial['magnetic_dot_induction']/MU0
    weight = geometry.weight[..., None]
    metrics = dict(current_relative_L2=relative(result.current_samples-exact['current'], exact['current'], weight),
                   induction_relative_L2=relative(rate-exact['induction'], exact['induction'], weight),
                   g_rate_relative_L2=relative(result.g_rate-initial['g_rate'], initial['g_rate'], mesh.volumes),
                   axis_g_rate_relative_L2=relative(result.g_rate[0]-initial['g_rate'][0], initial['g_rate'][0], np.ones(mesh.nz)),
                   psi_rate_max_relative=float(np.max(abs(result.psi_rate-initial['psi_rate']))/np.max(abs(initial['psi_rate']))) if np.max(abs(initial['psi_rate'])) > 0 else None,
                   local_energy_defect_relative_L2=relative(magnetic-total, local_reference, mesh.volumes),
                   local_energy_defect_rms_w_m3=float(np.sqrt(np.sum(mesh.volumes*(magnetic-total)**2)/np.sum(mesh.volumes))),
                   global_hall_work_relative=abs(result.magnetic_work_w)/result.absolute_work_scale_w,
                   magnetic_work_w=result.magnetic_work_w, work_scale_w=result.absolute_work_scale_w,
                   quadrature_reference_change=max(relative(initial[key]-control[key], control[key], np.ones_like(initial[key])) or 0. for key in ('g', 'g_rate', 'density')))
    arrays = dict(psi=initial['psi'], g=initial['g'], psi_rate=result.psi_rate, g_rate=result.g_rate,
                  reference_psi_rate=initial['psi_rate'], reference_g_rate=initial['g_rate'], magnetic_rate=magnetic,
                  total_energy_rate=total, electric_coefficients=result.electric_coefficients,
                  radial_faces=mesh.rf, axial_faces=mesh.zf)
    out = save(directory, case, dict(kind='spatial', elapsed_seconds=perf_counter()-start, metrics=metrics), arrays, source)
    print('Completed', case['name'], 'axis error', metrics['axis_g_rate_relative_L2'], flush=True)
    return out
