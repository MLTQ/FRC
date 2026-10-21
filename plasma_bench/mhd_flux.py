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


def rusanov(left, right, normal_field, axis, gamma, magnetic_dissipation=True, magnetic_energy=True, radius=None):
    """LLF face flux.

    magnetic_dissipation=False drops the tangential magnetic-energy jump from the energy
    dissipation; magnetic_energy=False also drops the magnetic part of the central energy flux,
    (B^2/mu0) v_n - (v.B) B_n/mu0. In either case the caller must carry the omitted magnetic
    energy transport as the Poynting flux of the constrained-transport electric field.

    2.5D face states (eight columns: ..., Omega = v_theta/r, g = B_theta/r) use the conservative
    angular-momentum and B_theta/r fluxes of `toroidal_flux.physical_flux` at the face `radius`;
    B^2 and v.B above then include the toroidal components.
    """
    left, right = left.copy(), right.copy()
    left[..., 4+axis] = normal_field
    right[..., 4+axis] = normal_field
    if left.shape[-1] == 8:
        return _toroidal_rusanov(left, right, axis, gamma, magnetic_dissipation, magnetic_energy, radius)
    ul, fl = physical_flux(left, axis, gamma)
    ur, fr = physical_flux(right, axis, gamma)
    speed = np.maximum(np.abs(left[..., 1+axis])+fast_speed(left, axis, gamma),
                       np.abs(right[..., 1+axis])+fast_speed(right, axis, gamma))
    jump = ur-ul
    if not (magnetic_dissipation and magnetic_energy):
        jump[..., 3] -= ((right[..., 4]**2+right[..., 5]**2)-(left[..., 4]**2+left[..., 5]**2))/(2*MU0)
    if not magnetic_energy:
        for q, f in ((left, fl), (right, fr)):
            vn, bn = q[..., 1+axis], q[..., 4+axis]
            f[..., 3] -= ((q[..., 4]**2+q[..., 5]**2)*vn-(q[..., 1]*q[..., 4]+q[..., 2]*q[..., 5])*bn)/MU0
    return .5*(fl+fr)-.5*speed[..., None]*jump, speed


def _toroidal_rusanov(left, right, axis, gamma, magnetic_dissipation, magnetic_energy, radius):
    from .toroidal_flux import physical_flux as toroidal_flux, fast_speed as toroidal_speed
    ul, fl = toroidal_flux(left, radius, axis, gamma)
    ur, fr = toroidal_flux(right, radius, axis, gamma)
    speed = np.maximum(np.abs(left[..., 1+axis])+toroidal_speed(left, radius, axis, gamma),
                       np.abs(right[..., 1+axis])+toroidal_speed(right, radius, axis, gamma))
    square = lambda q: q[..., 4]**2+q[..., 5]**2+(radius*q[..., 7])**2
    jump = ur-ul
    if not (magnetic_dissipation and magnetic_energy):
        jump[..., 3] -= (square(right)-square(left))/(2*MU0)
    if not magnetic_energy:
        for q, f in ((left, fl), (right, fr)):
            vn, bn = q[..., 1+axis], q[..., 4+axis]
            dot = q[..., 1]*q[..., 4]+q[..., 2]*q[..., 5]+radius**2*q[..., 6]*q[..., 7]
            f[..., 3] -= (square(q)*vn-dot*bn)/MU0
    return .5*(fl+fr)-.5*speed[..., None]*jump, speed


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
