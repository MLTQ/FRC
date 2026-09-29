"""Face energy flux for the curl-compatible Hall/resistive electric field."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .compatible_magnetic import radial_mass
from .toroidal_state import radius_squared


def radial_faces(mesh, value):
    """Centered magnetic trace in annular r-squared centroids; copy endpoints."""
    out = np.pad(value, ((1, 1), (0, 0)), mode='edge')[:-1].copy()
    s = radius_squared(mesh)[:, 0]
    fraction = ((mesh.rf[1:-1]**2-s[:-1])/(s[1:]-s[:-1]))[:, None]
    out[1:-1] = (1-fraction)*value[:-1]+fraction*value[1:]
    return out


def energy_flux(geometry, psi, g, coefficients):
    mesh = geometry.mesh
    er, et, ez = geometry.split_electric(coefficients)
    er = er.reshape(mesh.nr, mesh.nz)
    er = np.concatenate((er, er[:, :1]), axis=1)
    et = np.pad(et.reshape(mesh.nr-1, mesh.nz), ((1, 1), (0, 0)))
    et = np.concatenate((et, et[:, :1]), axis=1)
    ez = ez.reshape(-1, mesh.nz)
    if ez.shape[0] == mesh.nr:
        ez = np.pad(ez, ((0, 1), (0, 0)))
    elif ez.shape[0] != mesh.nr+1 or not hasattr(geometry, 'wall_btheta'):
        raise ValueError('Explicit magnetic wall trace required for free wall Ez')
    br, bz = mesh.face_fields(psi)
    bzf, gf = radial_faces(mesh, bz), radial_faces(mesh, g)
    if hasattr(geometry, 'wall_btheta'):
        gf[-1] = geometry.wall_btheta/mesh.radius
    product_z = (et[:, :-1]*bzf[:, :-1]+et[:, 1:]*bzf[:, 1:])/3
    product_z += (et[:, :-1]*bzf[:, 1:]+et[:, 1:]*bzf[:, :-1])/6
    sr = (product_z-ez*mesh.rf[:, None]*gf)/MU0
    gz = .5*(g+np.roll(g, 1, axis=1))
    gz = np.concatenate((gz, gz[:, :1]), axis=1)
    brz = .5*(br+np.roll(br, 1, axis=1))
    brz = np.concatenate((brz, brz[:, :1]), axis=1)
    ll, lr, rr = radial_mass(mesh.nr, mesh.dr)
    product_r = ll*et[:-1]*brz[:-1]+lr*(et[:-1]*brz[1:]+et[1:]*brz[:-1])+rr*et[1:]*brz[1:]
    sz = (radius_squared(mesh)*er*gz-product_r)/MU0
    return sr, sz


def energy_rate(geometry, psi, g, coefficients):
    sr, sz = energy_flux(geometry, psi, g, coefficients)
    mesh = geometry.mesh
    return -np.diff(mesh.rf[:, None]*sr, axis=0)/mesh.radial_measure[:, None]-np.diff(sz, axis=1)/mesh.dz
