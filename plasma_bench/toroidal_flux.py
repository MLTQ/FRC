"""Full three-component ideal MHD fluxes with conservative cylindrical weights."""
import numpy as np
from scipy.constants import mu_0 as MU0


def fast_speed(q, radius, axis, gamma):
    rho, p = q[..., 0], q[..., 3]
    b2 = (q[..., 4]**2+q[..., 5]**2+(radius*q[..., 7])**2)/(MU0*rho)
    bn2 = q[..., 4+axis]**2/(MU0*rho)
    a2 = gamma*p/rho
    return np.sqrt(.5*(a2+b2+np.sqrt(np.maximum((a2+b2)**2-4*a2*bn2, 0.))))


def physical_flux(q, radius, axis, gamma):
    rho, vr, vz, p, br, bz, omega, g = np.moveaxis(q, -1, 0)
    vt, bt = radius*omega, radius*g
    vn, bn = (vr, br) if axis == 0 else (vz, bz)
    benergy = (br**2+bz**2+bt**2)/(2*MU0)
    pt = p+benergy
    energy = p/(gamma-1)+.5*rho*(vr**2+vz**2+vt**2)+benergy
    angular = rho*radius**2*omega
    u = np.stack((rho, rho*vr, rho*vz, energy, angular, g), axis=-1)
    flux = np.stack((rho*vn, rho*vr*vn-br*bn/MU0, rho*vz*vn-bz*bn/MU0,
                     (energy+pt)*vn-(vr*br+vz*bz+vt*bt)*bn/MU0,
                     angular*vn-radius**2*g*bn/MU0, g*vn-omega*bn), axis=-1)
    flux[..., 1+axis] += pt
    return u, flux


def rusanov(left, right, normal_field, radius, axis, gamma):
    left, right = left.copy(), right.copy()
    left[..., 4+axis], right[..., 4+axis] = normal_field, normal_field
    ul, fl = physical_flux(left, radius, axis, gamma)
    ur, fr = physical_flux(right, radius, axis, gamma)
    speed = np.maximum(abs(left[..., 1+axis])+fast_speed(left, radius, axis, gamma),
                       abs(right[..., 1+axis])+fast_speed(right, radius, axis, gamma))
    return .5*(fl+fr)-.5*speed[..., None]*(ur-ul)
