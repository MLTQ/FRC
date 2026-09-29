"""Continuum-averaged initial states for nonlinear compressible Hall verification."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0
from plasma_bench.toroidal_state import ToroidalState, from_primitive
from .hall_work_snapshots import snapshot


def wave_state(mesh, reference, time=0., order=10):
    """Average true continuum conserved variables; do not fit recovered pressure."""
    x, w = leggauss(order)
    r = mesh.r[:, None, None, None]+.5*mesh.dr*x[None, None, :, None]
    z = mesh.z[None, :, None, None]+.5*mesh.dz*x[None, None, None, :]
    q = reference.fields(r, z, time)
    rho, pressure = reference.density+q['density'], reference.pressure+q['pressure']
    if np.any(rho <= 0) or np.any(pressure <= 0):
        raise ValueError('Wave amplitude violates the positive continuum reference domain')
    v, b = q['velocity'], q['magnetic']+np.array([0., 0., reference.bias])
    energy = pressure/(reference.gamma-1)+.5*rho*np.sum(v*v, axis=-1)+np.sum(b*b, axis=-1)/(2*MU0)
    fields = np.stack((rho, rho*v[..., 0], rho*v[..., 2], energy, rho*r*v[..., 1], b[..., 1]/r), axis=-1)
    weight = r*w[None, None, :, None]*w[None, None, None, :]*mesh.dr/(4*mesh.radial_measure[:, None, None, None])
    averages = np.sum(weight[..., None]*fields, axis=(2, 3))
    rc, zc = mesh.corners
    psi = .5*reference.bias*rc*rc+reference.fields(rc, zc, time)['psi']
    return ToroidalState(averages, psi, 'compatible')


def nonlinear_state(mesh, reference, ion_mass, pressure=16., gamma=5/3):
    data = snapshot(reference, mesh)
    return from_primitive(mesh, data['psi'], ion_mass*data['density'], pressure,
                          toroidal_over_radius=data['g'], gamma=gamma, magnetic_quadrature='compatible')
