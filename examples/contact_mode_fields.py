"""Independent physical reconstruction of the saved stationary growing mode."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.plasma_coordinates import PlasmaCoordinates
from plasma_bench.regular_current import RegularCurrent
from .regular_plasma_reference import initialize
from .current_transport_sampling import TransportSampler, bases
from .implicit_plasma_records import sha


def magnetic_curl(sample, magnetic):
    mesh = sample.mesh; count = (mesh.nr+1)*mesh.nz
    psi = magnetic[:count].reshape(mesh.nr+1, mesh.nz); g = magnetic[count:].reshape(mesh.nr, mesh.nz)
    dxx = bases(.5*mesh.rf**2, .5*sample.r**2, 5)[2]
    dzz = bases(mesh.zf, sample.z, 5, True)[2]
    return np.stack((-sample.R*(sample.ex@g@sample.dez.T),
        -(sample.nx@psi@dzz.T)/sample.R-sample.R*(dxx@psi@sample.nz.T),
        (2*sample.ex+sample.R**2*sample.dex)@g@sample.ez.T), axis=-1)


def reconstruct(order=16):
    directory = Path('outputs/stationary_contact_linear_probe')
    record = json.loads((directory/'case_02.json').read_text()); path = directory/record['arrays_file']
    if sha(path) != record['arrays_sha256']: raise ValueError('Original frozen mode required')
    raw = dict(np.load(path))
    ref, d, _, system, field = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    u = raw['state']; c = PlasmaCoordinates(system, u, field)
    if record['leading_imaginary'] != 0 or np.max(abs(raw['leading_vector'].imag)) != 0:
        raise ValueError('This diagnostic expects the saved real growing mode')
    du, db = c.unpack(raw['leading_vector'].real)
    sample = TransportSampler(d, system.geometry.mesh, order)
    n, _, _, dn, velocity = sample.fluid(u, du)
    entropy0 = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
    entropy = np.stack([(sample.scalar(du[..., k])-entropy0[..., k-4]*dn)/n for k in (4, 5)], axis=-1)
    pressure0 = E*n[..., None]*np.exp(entropy0)*(n[..., None]/d.thermal.n0)**(d.thermal.gammas-1)
    model = RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)
    current_coefficients = model.current(db); current, gradient = sample.electric(current_coefficients, True)
    magnetic = sample.magnetic(db); mesh = model.mesh; count = (mesh.nr+1)*mesh.nz
    g = db[count:].reshape(mesh.nr, mesh.nz); curl_b = magnetic_curl(sample, db)
    div_j = np.trace(gradient, axis1=-2, axis2=-1)
    # Exact one-sided outer trace: no offset from the physical boundary.
    edge_wall = bases(.5*mesh.rf**2, np.array([.5*mesh.radius**2]), 5)[3]
    jr = current_coefficients[:mesh.nr*mesh.nz].reshape(mesh.nr, mesh.nz)
    wall_jr = mesh.radius*(edge_wall@jr@sample.nz.T).ravel()
    wall_btheta = mesh.radius*(edge_wall@g@sample.ez.T).ravel()
    wall_vr = mesh.radius*((np.fft.fft(du[-1, -1, :, 1]/u[-1, -1, :, 0])/d.nz)@sample.fourier.T).real
    radial_weight = sample.weight.sum(axis=1)/d.length
    axial_weight = sample.weight.sum(axis=0)/(np.pi*d.radius**2)
    integral = lambda value: float(np.sum(sample.weight*value))
    j2 = integral(np.sum(current*current, axis=-1)); volume = float(sample.weight.sum())
    target = ref.fields(sample.Z)
    target_pressure = np.stack(np.broadcast_arrays(target['electron_pressure'], target['ion_pressure']), axis=-1)
    result = dict(order=order, eigenvalue=record['maximum_real'], input_arrays_sha256=sha(path),
        strong_ampere_relative=float(np.sqrt(integral(np.sum((current-curl_b/MU0)**2, axis=-1))/j2)),
        strong_ampere_components=[float(np.sqrt(integral((current[..., k]-curl_b[..., k]/MU0)**2)/integral(current[..., k]**2))) for k in range(3)],
        strong_charge_relative=float(np.sqrt(integral(div_j*div_j)/integral(np.sum(gradient*gradient, axis=(-2, -1))))),
        wall_current_rms_over_volume_rms=float(np.sqrt(np.sum(axial_weight*wall_jr**2)/d.length/(j2/volume))),
        wall_velocity_rms_over_radial_volume_rms=float(np.sqrt(np.sum(axial_weight*wall_vr**2)/d.length/(integral(velocity[..., 0]**2)/volume))),
        wall_btheta_rms_over_volume_rms=float(np.sqrt(np.sum(axial_weight*wall_btheta**2)/d.length/(integral(magnetic[..., 1]**2)/volume))),
        background_pressure_relative=[float(np.sqrt(integral((pressure0[..., k]-target_pressure[..., k])**2)/integral(np.broadcast_to(target_pressure[..., k]**2, n.shape)))) for k in range(2)],
        background_total_pressure_max_relative=float(np.max(abs(pressure0.sum(axis=-1)-ref.total_pressure))/ref.total_pressure))
    physical = np.concatenate((dn[..., None], velocity, entropy, magnetic), axis=-1)
    result['outer_field_cell_norm_fractions'] = (np.sum(sample.weight[..., None]*(sample.R > mesh.rf[-2])[..., None]*physical**2, axis=(0, 1))/np.sum(sample.weight[..., None]*physical**2, axis=(0, 1))).tolist()
    arrays = dict(radius=sample.r, axial=sample.z, radial_weight=radial_weight, axial_weight=axial_weight,
        physical_mode=physical, current=current, curl_magnetic=curl_b, current_gradient=gradient,
        charge_divergence=div_j, wall_current=wall_jr, wall_btheta=wall_btheta, wall_radial_velocity=wall_vr, density=n, pressure=pressure0,
        fluid_change=du, magnetic_change=db, current_coefficients=current_coefficients)
    return result, arrays
