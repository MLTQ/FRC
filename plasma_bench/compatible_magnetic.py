"""Positive magnetic energy of a curl-compatible cylindrical flux interpolant."""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0


@lru_cache(maxsize=16)
def radial_mass(nr, dr):
    """Annular means of radial shape-function products; fixed accurate quadrature."""
    faces = np.arange(nr+1)*dr
    x, w = leggauss(16)
    r = (faces[:-1, None]+faces[1:, None])/2+dr*x/2
    measure = np.diff(faces**2)/2
    fraction = (r*r-faces[:-1, None]**2)/(2*measure[:, None])
    left = faces[:-1, None]/r*(1-fraction)
    right = faces[1:, None]/r*fraction
    weight = r*w*dr/(2*measure[:, None])
    values = tuple(np.sum(weight*v, axis=1)[:, None] for v in (left*left, left*right, right*right))
    for value in values:
        value.setflags(write=False)
    return values


def energy_density(mesh, psi):
    """Poloidal energy of psi linear in x=r²/2 and z within each cell."""
    br, bz = mesh.face_fields(psi)
    ll, lr, rr = radial_mass(mesh.nr, mesh.dr)
    radial = ll*br[:-1]**2+2*lr*br[:-1]*br[1:]+rr*br[1:]**2
    axial = (bz[:, :-1]**2+bz[:, :-1]*bz[:, 1:]+bz[:, 1:]**2)/3
    return (radial+axial)/(2*MU0)


def energy_rate(mesh, psi, psi_rate):
    """Exact directional derivative of this quadratic cell energy."""
    br, bz = mesh.face_fields(psi)
    dr, dz = mesh.face_fields(psi_rate)
    ll, lr, rr = radial_mass(mesh.nr, mesh.dr)
    radial = (ll*br[:-1]+lr*br[1:])*dr[:-1]+(lr*br[:-1]+rr*br[1:])*dr[1:]
    axial = ((bz[:, :-1]+.5*bz[:, 1:])*dz[:, :-1]+(.5*bz[:, :-1]+bz[:, 1:])*dz[:, 1:])/3
    return (radial+axial)/MU0


def gradient(mesh, psi):
    """Full periodic corner gradient of the integrated poloidal energy."""
    br, bz = mesh.face_fields(psi)
    ll, lr, rr = radial_mass(mesh.nr, mesh.dr)
    ar = np.zeros_like(br)
    ar[:-1] += mesh.volumes*(ll*br[:-1]+lr*br[1:])/MU0
    ar[1:] += mesh.volumes*(lr*br[:-1]+rr*br[1:])/MU0
    az = np.zeros_like(bz)
    az[:, :-1] += mesh.volumes*(bz[:, :-1]+.5*bz[:, 1:])/(3*MU0)
    az[:, 1:] += mesh.volumes*(.5*bz[:, :-1]+bz[:, 1:])/(3*MU0)
    radial = np.divide(ar, mesh.rf[:, None]*mesh.dz, out=np.zeros_like(ar), where=mesh.rf[:, None] > 0)
    axial = az/mesh.radial_measure[:, None]
    result = np.zeros_like(psi)
    result[:, :-1] += radial; result[:, 1:] -= radial
    result[:-1] -= axial; result[1:] += axial
    result[:, 0] += result[:, -1]
    result[:, -1] = 0.  # The periodic duplicate is not a second independent DOF.
    return result
