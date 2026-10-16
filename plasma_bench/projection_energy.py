"""Energy transfer paired with a Galerkin momentum projection."""
import numpy as np


def projected_frame_energy(energy, left_momentum, right_momentum, left_velocity,
                           right_velocity, kinetic):
    """Return thermal energy and the explicit local total-energy redistribution.

    Let W_c=integral(p_source dot v_target) and K_c be target kinetic energy.
    U_c=E_source,c-W_c+K_c, and delta E_c=2 K_c-W_c.
    The summed redistribution is v^T(Mv-b), hence zero for the exact projection.
    Nothing is clipped, floored, or adjusted using a global energy residual.
    """
    work = np.sum(left_momentum*left_velocity+right_momentum*right_velocity, axis=1)
    thermal = energy-work+kinetic
    redistribution = 2*kinetic-work
    if not np.isfinite(thermal).all() or not np.isfinite(redistribution).all():
        raise ValueError('Unrepresentable projected-frame energy')
    return thermal, redistribution
