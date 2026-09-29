"""Independent conservative averages of the analytic nonlinear Hall snapshot."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.special import j0, j1


def snapshot(reference, mesh, order=10):
    x, w = leggauss(order)
    r = mesh.r[:, None, None, None]+.5*mesh.dr*x[None, None, :, None]
    z = mesh.z[None, :, None, None]+.5*mesh.dz*x[None, None, None, :]
    data = reference.fields(r, z)
    weights = w[None, None, :, None]*w[None, None, None, :]
    factor = mesh.dr/(4*mesh.radial_measure[:, None])
    # g is the annular average of Btheta/r; its r in the volume measure cancels.
    g = np.sum(weights*data['magnetic'][..., 1], axis=(-2, -1))*factor
    g_rate = np.sum(weights*data['induction'][..., 1], axis=(-2, -1))*factor
    density = np.sum(weights*r*data['number_density'], axis=(-2, -1))*factor
    magnetic_rate = np.sum(data['magnetic']*data['induction'], axis=-1)
    rate_mean = np.sum(weights*r*magnetic_rate, axis=(-2, -1))*factor
    corner = reference.fields(*mesh.corners)
    return dict(psi=corner['psi'], g=g, g_rate=g_rate, psi_rate=corner['psi_rate'],
                density=density, magnetic_dot_induction=rate_mean)


def wave_hall_rate(reference, r, z):
    """Independently differentiated linear Hall part, without fluid induction."""
    r, z = np.broadcast_arrays(r, z)
    bp, bt = reference.eigenpair[1][4:]
    a = np.sqrt(MU0*reference.density)*bp/reference.total_k
    t = np.sqrt(MU0*reference.density)*bt
    da = -1j*reference.k*reference.bias*t/(MU0*reference.number_density*E)
    dt = -1j*reference.k*reference.bias*reference.total_k**2*a/(MU0*reference.number_density*E)
    phase = np.exp(1j*reference.k*z)
    f1, f0 = j1(reference.q*r), j0(reference.q*r)
    magnetic = np.stack((-1j*reference.k*da*f1*phase, dt*f1*phase, reference.q*da*f0*phase), axis=-1).real
    return magnetic, (r*da*f1*phase).real


def wave_snapshot(reference, mesh, order=10):
    x, w = leggauss(order)
    r = mesh.r[:, None, None, None]+.5*mesh.dr*x[None, None, :, None]
    z = mesh.z[None, :, None, None]+.5*mesh.dz*x[None, None, None, :]
    data = reference.fields(r, z)
    rate, _ = wave_hall_rate(reference, r, z)
    weights = w[None, None, :, None]*w[None, None, None, :]
    factor = mesh.dr/(4*mesh.radial_measure[:, None])
    g = np.sum(weights*data['magnetic'][..., 1], axis=(-2, -1))*factor
    g_rate = np.sum(weights*rate[..., 1], axis=(-2, -1))*factor
    magnetic = data['magnetic']+np.array([0., 0., reference.bias])
    rate_mean = np.sum(weights*r*np.sum(magnetic*rate, axis=-1), axis=(-2, -1))*factor
    rc, zc = mesh.corners
    _, psi_rate = wave_hall_rate(reference, rc, zc)
    return dict(psi=.5*reference.bias*rc*rc+reference.fields(rc, zc)['psi'], g=g, g_rate=g_rate,
                psi_rate=psi_rate, density=np.full_like(g, reference.number_density), magnetic_dot_induction=rate_mean)
