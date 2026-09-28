"""Rusanov fluid fluxes using a single constrained normal magnetic field."""
import numpy as np
from scipy.constants import mu_0 as MU0


def fast_speed(q, axis, gamma):
    rho, p = q[..., 0], q[..., 3]
    b2 = (q[..., 4]**2+q[..., 5]**2)/(MU0*rho)
    bn2 = q[..., 4+axis]**2/(MU0*rho)
    a2 = gamma*p/rho
    discriminant = np.maximum((a2+b2)**2-4*a2*bn2, 0)
    return np.sqrt(.5*(a2+b2+np.sqrt(discriminant)))


def physical_flux(q, axis, gamma):
    rho, vr, vz, p, br, bz = np.moveaxis(q, -1, 0)
    vn = vr if axis == 0 else vz
    bn = br if axis == 0 else bz
    pt = p+(br*br+bz*bz)/(2*MU0)
    energy = p/(gamma-1)+.5*rho*(vr*vr+vz*vz)+(br*br+bz*bz)/(2*MU0)
    u = np.stack((rho, rho*vr, rho*vz, energy), axis=-1)
    flux = np.stack((rho*vn, rho*vr*vn-br*bn/MU0,
                     rho*vz*vn-bz*bn/MU0, (energy+pt)*vn-(vr*br+vz*bz)*bn/MU0), axis=-1)
    flux[..., 1+axis] += pt
    return u, flux


def rusanov(left, right, normal_field, axis, gamma):
    left, right = left.copy(), right.copy()
    left[..., 4+axis] = normal_field
    right[..., 4+axis] = normal_field
    ul, fl = physical_flux(left, axis, gamma)
    ur, fr = physical_flux(right, axis, gamma)
    speed = np.maximum(np.abs(left[..., 1+axis])+fast_speed(left, axis, gamma),
                       np.abs(right[..., 1+axis])+fast_speed(right, axis, gamma))
    return .5*(fl+fr)-.5*speed[..., None]*(ur-ul), speed


def ghost_cells(q):
    """Axis reflection, outer radial slip wall, and periodic axial ghost cells."""
    out = np.pad(q, ((1, 1), (1, 1), (0, 0)), mode="edge")
    out[0, 1:-1, 1] *= -1
    out[0, 1:-1, 4] *= -1
    out[-1, 1:-1, 1] *= -1
    out[-1, 1:-1, 4] *= -1
    out[:, 0] = out[:, -2]
    out[:, -1] = out[:, 1]
    return out
