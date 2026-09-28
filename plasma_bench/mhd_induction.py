"""Upwind corner electric fields and explicitly resistive flux evolution."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_flux import fast_speed
from .mhd_reconstruction import staggered_corner_fields


def ideal_emf(mesh, ghost, br_faces, bz_faces, gamma, reconstruction=None):
    if reconstruction is None:
        sw, se = ghost[:-1, :-1], ghost[1:, :-1]
        nw, ne = ghost[:-1, 1:], ghost[1:, 1:]
        bz_pad = np.pad(bz_faces, ((1, 1), (0, 0)), mode="edge")
        bz_left, bz_right = bz_pad[:-1], bz_pad[1:]
        br_pad = np.concatenate((br_faces[:, -1:], br_faces, br_faces[:, :1]), axis=1)
        br_down, br_up = br_pad[:, :-1], br_pad[:, 1:]
    else:
        sw, se, nw, ne = reconstruction.corners()
        bz_left, bz_right, br_down, br_up = staggered_corner_fields(mesh, br_faces, bz_faces)
    ar = np.maximum.reduce([np.abs(q[..., 1])+fast_speed(q, 0, gamma) for q in (sw, se, nw, ne)])
    az = np.maximum.reduce([np.abs(q[..., 2])+fast_speed(q, 1, gamma) for q in (sw, se, nw, ne)])
    electric = .25*((sw[..., 1]+nw[..., 1])*bz_left+(se[..., 1]+ne[..., 1])*bz_right
                     -(sw[..., 2]+se[..., 2])*br_down-(nw[..., 2]+ne[..., 2])*br_up)
    electric += -.5*ar*(bz_right-bz_left)+.5*az*(br_up-br_down)
    electric[[0, -1]] = 0  # Axis regularity / fixed perfectly conducting wall.
    electric[:, -1] = electric[:, 0]
    return electric


def resistive_emf(mesh, psi, resistivity):
    """eta Jtheta at corners; fixed boundary flux and periodic axial derivative."""
    out = np.zeros_like(psi)
    if resistivity == 0:
        return out
    unique = psi[:, :-1]
    rr = mesh.rf[1:-1, None]
    radial = (unique[2:]-2*unique[1:-1]+unique[:-2])/mesh.dr**2
    radial -= (unique[2:]-unique[:-2])/(2*rr*mesh.dr)
    axial = (np.roll(unique[1:-1], -1, axis=1)-2*unique[1:-1]+np.roll(unique[1:-1], 1, axis=1))/mesh.dz**2
    out[1:-1, :-1] = -resistivity*(radial+axial)/(MU0*rr)
    out[:, -1] = out[:, 0]
    return out


def resistive_energy_flux(mesh, q, electric):
    # Average tangential B to the fluid face. Shared electric field ensures no
    # unreported energy source: only divergence of the resistive Poynting flux.
    bz_r = np.pad(q[..., 5], ((1, 1), (0, 0)), mode="edge")
    bz_r = .5*(bz_r[:-1]+bz_r[1:])
    br_z = np.concatenate((q[:, -1:, 4], q[..., 4], q[:, :1, 4]), axis=1)
    br_z = .5*(br_z[:, :-1]+br_z[:, 1:])
    return .5*(electric[:, :-1]+electric[:, 1:])*bz_r/MU0, -.5*(electric[:-1]+electric[1:])*br_z/MU0


def ohmic_power(mesh, electric, resistivity):
    """Corner quadrature of eta J²; an independent physical-heating estimate."""
    if resistivity == 0:
        return 0.0
    # Interior radial corner dual volumes, and unique periodic z corners.
    return float(np.sum(electric[1:-1, :-1]**2/resistivity *
                        (2*np.pi*mesh.rf[1:-1, None]*mesh.dr*mesh.dz)))
