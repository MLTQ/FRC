"""Explicit face transport audit for the experimental Hall electric projection."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0
from .toroidal_state import radius_squared


def energy_flux(geometry, psi, g, electric_coefficients):
    """Return S_r/S_z face averages; no magnetic or thermal energy correction."""
    mesh = geometry.mesh
    er, et, ez = geometry.split_electric(electric_coefficients)
    er = er.reshape(mesh.nr, mesh.nz)
    er = np.concatenate((er, er[:, :1]), axis=1)
    et = np.pad(et.reshape(mesh.nr-1, mesh.nz), ((1, 1), (0, 0)))
    et = np.concatenate((et, et[:, :1]), axis=1)
    ez = np.pad(ez.reshape(mesh.nr, mesh.nz), ((0, 1), (0, 0)))
    br, bz = mesh.face_fields(psi)
    def radial_face(values):
        out = np.pad(values, ((1, 1), (0, 0)), mode='edge')[:-1].copy()
        s = radius_squared(mesh)[:, 0]
        weight = ((mesh.rf[1:-1]**2-s[:-1])/(s[1:]-s[:-1]))[:, None]
        out[1:-1] = (1-weight)*values[:-1]+weight*values[1:]
        return out
    bzf, gf = radial_face(bz), radial_face(g)
    product_z = (et[:, :-1]*bzf[:, :-1]+et[:, 1:]*bzf[:, 1:])/3
    product_z += (et[:, :-1]*bzf[:, 1:]+et[:, 1:]*bzf[:, :-1])/6
    sr = (product_z-ez*mesh.rf[:, None]*gf)/MU0
    gz = .5*(g+np.roll(g, 1, axis=1))
    gz = np.concatenate((gz, gz[:, :1]), axis=1)
    brz = .5*(br+np.roll(br, 1, axis=1))
    brz = np.concatenate((brz, brz[:, :1]), axis=1)
    x, w = leggauss(3)
    r = mesh.r[:, None, None]+.5*mesh.dr*x[None, None, :]
    fraction = .5*(x+1)[None, None, :]
    theta = (1-fraction)*et[:-1, :, None]+fraction*et[1:, :, None]
    radial = (1-fraction)*brz[:-1, :, None]+fraction*brz[1:, :, None]
    product_r = np.sum(r*theta*radial*w, axis=-1)*mesh.dr/(2*mesh.radial_measure[:, None])
    sz = (radius_squared(mesh)*er*gz-product_r)/MU0
    return sr, sz


def local_energy_audit(geometry, psi, g, result):
    """Magnetic rate and total-energy flux rate are independent, not repaired."""
    mesh = geometry.mesh
    sr, sz = energy_flux(geometry, psi, g, result.electric_coefficients)
    total_rate = -np.diff(mesh.rf[:, None]*sr, axis=0)/mesh.radial_measure[:, None]-np.diff(sz, axis=1)/mesh.dz
    br, bz = mesh.cell_fields(psi)
    dbr, dbz = mesh.cell_fields(result.psi_rate)
    magnetic_rate = (br*dbr+bz*dbz+radius_squared(mesh)*g*result.g_rate)/MU0
    return magnetic_rate, total_rate, sr, sz
