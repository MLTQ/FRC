"""Polarization-preserving electric-dipole emission for a J=0-1-0 line."""
import numpy as np
from .resonance_transport import isotropic


def dipole_emission(rng, dipoles):
    """Sample 3[1-(d.n)^2]/(8*pi); return direction and transverse polarization.

    Each real unit dipole represents one linearly polarized pure state. A
    statistical mixture carries unpolarized or partially linearly polarized
    light without resetting that mixture at a subsequent absorption.
    """
    dipoles = np.asarray(dipoles, float)
    if (dipoles.ndim != 2 or dipoles.shape[1] != 3 or not np.isfinite(dipoles).all()
            or np.any(abs(np.sum(dipoles*dipoles, axis=1)-1) > 2e-12)):
        raise ValueError('Finite real unit dipoles required')
    outgoing = np.empty_like(dipoles); polarization = np.empty_like(dipoles)
    remaining = np.arange(len(dipoles))
    while len(remaining):
        direction = isotropic(rng, len(remaining))
        projected = dipoles[remaining]-direction*np.sum(dipoles[remaining]*direction, axis=1)[:, None]
        norm2 = np.sum(projected*projected, axis=1)
        accept = rng.random(len(remaining)) < norm2
        ids = remaining[accept]
        outgoing[ids] = direction[accept]
        polarization[ids] = projected[accept]/np.sqrt(norm2[accept])[:, None]
        remaining = remaining[~accept]
    return outgoing, polarization


def observer_stokes(dipoles):
    """I,Q,U per steradian along +x, with positive Q along laboratory +z."""
    y, z = np.asarray(dipoles)[:, 1:].T
    return 3/(8*np.pi)*np.column_stack((y*y+z*z, z*z-y*y, 2*y*z))
