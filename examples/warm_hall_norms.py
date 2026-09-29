"""Physical species errors for warm Hall waves and continuous Q1 comparisons."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.toroidal_state import primitive
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction
from .compatible_hall_evolution import wave_errors as fluid_wave_errors
from .warm_hall_spatial import relative


def wave_errors(model, reference):
    metrics, exact_state = fluid_wave_errors(model, reference, model.time)
    geo, mesh = model.geometry, model.mesh
    pe = geo.scalar_samples(model.thermal.pressure(model.state, model.electron_entropy))
    ne = geo.scalar_samples(model.thermal.density(model.state))
    rec = ToroidalReconstruction(mesh, primitive(mesh, model.state, reference.gamma), order=model.parameters.spatial_order)
    total, rho = np.empty_like(geo.r), np.empty_like(geo.r)
    for i in range(geo.order):
        r = geo.r[:, 0, i, 0]
        for j in range(geo.order):
            offset = (geo.z[0, 0, 0, j]-mesh.z[0])/mesh.dz
            q = rec.primitive(rec.raw(r, offset), rec.auxiliary(r, offset), r[:, None])
            rho[:, :, i, j], total[:, :, i, j] = q[..., 0], q[..., 3]
    pi = total-pe
    f = reference.fields(geo.r, geo.z, model.time)
    target_pe = reference.electron_pressure+f['electron_pressure']
    target_pi = reference.pressure-reference.electron_pressure+f['ion_pressure']
    target_n = reference.number_density+f['density']/reference.ion_mass
    te, ti = pe/(ne*E), pi/(rho/reference.ion_mass*E)
    true_te, true_ti = target_pe/(target_n*E), target_pi/(target_n*E)
    for name, value, exact, variation in (
            ('electron_pressure', pe, target_pe, f['electron_pressure']),
            ('ion_pressure', pi, target_pi, f['ion_pressure']),
            ('electron_temperature', te, true_te, true_te-reference.electron_pressure/(reference.number_density*E)),
            ('ion_temperature', ti, true_ti, true_ti-(reference.pressure-reference.electron_pressure)/(reference.number_density*E))):
        metrics[name+'_relative_perturbation_L2'] = relative(value-exact, variation, geo.weight)
        metrics[name+'_rms_error'] = float(np.sqrt(np.sum(geo.weight*(value-exact)**2)/np.sum(geo.weight)))
    nodal_r, nodal_z = model.thermal.mapping.geometry.r, model.thermal.mapping.geometry.z
    c = reference.fields(nodal_r, nodal_z)['electron_pressure']
    s = reference.fields(nodal_r, nodal_z, np.pi/(2*abs(reference.omega)))['electron_pressure']
    actual = model.thermal.pressure(model.state, model.electron_entropy)-reference.electron_pressure
    weight = model.thermal.mapping.volumes
    gram = np.array([[np.sum(weight*a*b) for b in (c, s)] for a in (c, s)])
    coefficients = np.linalg.solve(gram, [np.sum(weight*actual*a) for a in (c, s)])
    metrics['electron_mode_amplitude'] = float(np.linalg.norm(coefficients))
    phase = np.arctan2(coefficients[1], coefficients[0])-abs(reference.omega)*model.time
    metrics['electron_mode_phase_error_rad'] = float(np.angle(np.exp(1j*phase)))
    return metrics, exact_state


def sample_q1(mesh, value, x, z):
    x, z = np.broadcast_arrays(x, z)
    edges = .5*mesh.rf**2
    i = np.clip(np.searchsorted(edges, x, side='right')-1, 0, mesh.nr-1)
    position = ((z+mesh.half_length)/mesh.dz) % mesh.nz
    j = np.floor(position).astype(int); tz = position-j
    tx = (x-edges[i])/(edges[i+1]-edges[i])
    return (1-tx)*((1-tz)*value[i, j]+tz*value[i, (j+1)%mesh.nz])+tx*((1-tz)*value[i+1, j]+tz*value[i+1, (j+1)%mesh.nz])


def pressure_difference(coarse_mesh, coarse, fine_mesh, fine, initial_fine):
    """Exact polynomial norm on nested fine cells, normalized by fine change."""
    nodes, weights = np.polynomial.legendre.leggauss(3)
    lo, hi = .5*fine_mesh.rf[:-1]**2, .5*fine_mesh.rf[1:]**2
    numerator = denominator = 0.
    for a, wa in zip(nodes, weights):
        x = .5*(lo+hi+a*(hi-lo))[:, None]
        for b, wb in zip(nodes, weights):
            z = fine_mesh.z[None]+.5*fine_mesh.dz*b
            ca = sample_q1(coarse_mesh, coarse, x, z)
            fi = sample_q1(fine_mesh, fine, x, z)
            initial = sample_q1(fine_mesh, initial_fine, x, z)
            weight = .25*wa*wb*fine_mesh.volumes
            numerator += np.sum(weight*(ca-fi)**2)
            denominator += np.sum(weight*(fi-initial)**2)
    return dict(rms_difference_pa=float(np.sqrt(numerator/np.sum(fine_mesh.volumes))),
                relative_to_fine_pressure_change=float(np.sqrt(numerator/denominator)) if denominator else None)
