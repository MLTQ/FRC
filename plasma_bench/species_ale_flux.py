"""Conservative species flux through a moving face in fixed velocity coordinates."""
import numpy as np
from .species_euler_flux import hllc


def moving_flux(left, right, mass_ratio, frame_normal, face_speed):
    """Dimensionless axial ALE flux; face_speed is measured in units of v0."""
    speed = np.asarray(face_speed)
    if speed.shape != left.shape[:1] or not np.isfinite(speed).all():
        raise ValueError('One finite moving speed per face required')
    l = left.copy(); r = right.copy()
    l[..., 4] -= speed[:, None]; r[..., 4] -= speed[:, None]
    flux = hllc(l, r, mass_ratio, frame_normal)
    # Transform face-rest-frame conservative fluxes back to the one fixed
    # velocity frame. In particular, retain pressure work on the moving face.
    result = flux.copy(); a = speed[:, None]
    result[..., 3] += mass_ratio*a*flux[..., 0]
    result[..., 4] += a*flux[..., 3]+.5*mass_ratio*a*a*flux[..., 0]
    return result
