"""Physical local pressure work and zero-Btheta boundary measurements."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.special import j1, jn_zeros
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_mass import DualMass
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_energy import energy_flux, energy_rate
from plasma_bench.hall_trace import TraceHallCurl
from plasma_bench.hall_midpoint import HallMidpoint
from .electron_reference import ElectronReference
from .pressure_work_reference import fields, cell_average, face_integrals


def relative(error, reference, weight):
    denominator = np.sum(weight*reference**2)
    return float(np.sqrt(np.sum(weight*error**2)/denominator)) if denominator else None


def pressure_spatial(case):
    ref = ElectronReference(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    mapping = DualMass(mesh); geo = mapping.geometry
    f = fields(ref, geo.r, geo.z, barotropic=case['barotropic'])
    model = ElectronPressure(mesh, f['density'], gamma=ref.gamma)
    u, g = model.from_pressure(f['pressure']), ref.cell_average(mesh, 'g')
    rate, flux = energy_rate(model, g, u), energy_flux(model, g, u)
    dg, du, _, _ = model.rhs(g, u)
    magnetic = geo.magnetic_mass*g*dg/mesh.volumes
    thermal = mapping.primal_average(model.energy_weight*du)
    exact = cell_average(ref, mesh, 'energy_rate', barotropic=case['barotropic'])
    exact_check = cell_average(ref, mesh, 'energy_rate', order=18, barotropic=case['barotropic'])
    work = cell_average(ref, mesh, 'pressure_work', barotropic=case['barotropic'])
    exact_flux = face_integrals(ref, mesh, barotropic=case['barotropic'])
    metrics = dict(energy_rate_relative_L2=relative(rate-exact, exact, mesh.volumes),
                   energy_error_relative_to_pressure_work=relative(rate-exact, work, mesh.volumes),
                   local_ion_defect_relative_to_pressure_work=relative(rate-magnetic-thermal, work, mesh.volumes),
                   local_ion_defect_rms_w_m3=float(np.sqrt(np.sum(mesh.volumes*(rate-magnetic-thermal)**2)/np.sum(mesh.volumes))),
                   axis_energy_error_relative_to_pressure_work=relative(rate[0]-exact[0], work[0], np.ones(mesh.nz)),
                   radial_face_flux_relative_L2=relative(flux[0]-exact_flux[0], exact_flux[0], np.ones_like(flux[0])),
                   axial_face_flux_relative_L2=relative(flux[1]-exact_flux[1], exact_flux[1], np.ones_like(flux[1])),
                   reference_quadrature_change=relative(exact-exact_check, work, mesh.volumes),
                   global_energy_rate_w=float(np.sum(mesh.volumes*rate)))
    arrays = dict(g=g, entropy=u, number_density=model.density, pressure=f['pressure'], energy_rate=rate,
                  magnetic_rate=magnetic, electron_rate=thermal, reference_energy_rate=exact, reference_pressure_work=work,
                  radial_flux=flux[0], axial_flux=flux[1], reference_radial_flux=exact_flux[0], reference_axial_flux=exact_flux[1],
                  primal_volumes=mesh.volumes, dual_volumes=geo.volumes, radial_faces=mesh.rf, axial_faces=mesh.zf)
    return dict(metrics=metrics), arrays


def trace_diffusion(case):
    ref = ElectronReference(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    geo = TraceHallCurl(mesh, case['quadrature'])
    q = jn_zeros(1, 1)[0]/mesh.radius
    nodes, weights = np.polynomial.legendre.leggauss(20)
    r = mesh.r[:, None]+.5*mesh.dr*nodes
    values = .001*np.sum(weights*j1(q*r), axis=1)*mesh.dr/(2*mesh.radial_measure)
    g = np.broadcast_to(values[:, None], (mesh.nr, mesh.nz)).copy()
    psi = np.zeros((mesh.nr+1, mesh.nz+1)); eta = case['resistivity']
    current = geo.current(psi, g)
    coefficients = geo.project_electric(eta*current)
    rate = geo.unpack(geo.faraday @ coefficients)[1]
    expected_rate = -eta/MU0*q*q*g
    end = case['decay_exponent']*MU0/(eta*q*q)
    initial_g, initial_psi = g.copy(), psi.copy()
    stepper = HallMidpoint(geo, tolerance=case['tolerance'], max_iterations=case['hall_iterations'])
    steps = case['steps'] or 1; heat = 0.; total = np.zeros_like(g)
    for _ in range(steps):
        result = stepper.step(psi, g, None, end/steps, resistivity=eta)
        psi, g = result.psi, result.g
        heat += result.physical_heat_j; total += result.total_energy_change
    exact = initial_g*np.exp(-case['decay_exponent'])
    s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
    initial_h = .5*np.sum(mesh.volumes*s*initial_g**2)/MU0
    final_h = np.sum(mesh.volumes*geo.energy_density(psi, g))
    metrics = dict(instantaneous_g_rate_relative_L2=relative(rate-expected_rate, expected_rate, mesh.volumes),
                   axis_instantaneous_g_rate_relative_signed_error=float(rate[0, 0]/expected_rate[0, 0]-1),
                   evolved_g_relative_L2=relative(g-exact, exact, mesh.volumes),
                   evolved_axis_relative_signed_error=float(g[0, 0]/exact[0, 0]-1),
                   evolved_magnetic_coefficient_relative_L2=relative(g-exact, exact, mesh.volumes*s),
                   energy_heat_relative_residual=float((final_h-initial_h+heat)/initial_h),
                   conservative_total_energy_change_j=float(np.sum(mesh.volumes*total)))
    return dict(metrics=metrics, end_time_s=end, accepted_steps=steps, physical_heat_j=heat), dict(initial_g=initial_g, final_g=g,
                initial_psi=initial_psi, final_psi=psi, reference_g=exact, g_rate=rate, reference_g_rate=expected_rate,
                total_energy_change=total, radial_faces=mesh.rf, axial_faces=mesh.zf, primal_volumes=mesh.volumes)
