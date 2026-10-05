"""Atomic classical RK4 for experiments with the bounded coupled equations."""
import numpy as np
from .harmonic_magnetic import HarmonicMagnetic


def field_from_packed(geometry, packed, bias):
    psi, g = geometry.unpack(packed)
    return HarmonicMagnetic(geometry.mesh, psi, g, bias)


def rk4(system, state, field, dt):
    """Return new fluid/field plus stage diagnostics; no stability step is inferred."""
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('Positive finite coupled time step required')
    state = np.asarray(state, float)
    magnetic = field.packed_increment(system.geometry)
    stages = []

    def rate(u, b):
        result = system.derivative(u, field_from_packed(system.geometry, b, field.axial_bias))
        stages.append(result)
        return result.fluid, result.magnetic

    k1, b1 = rate(state, magnetic)
    k2, b2 = rate(state+.5*dt*k1, magnetic+.5*dt*b1)
    k3, b3 = rate(state+.5*dt*k2, magnetic+.5*dt*b2)
    k4, b4 = rate(state+dt*k3, magnetic+dt*b3)
    updated = state+dt*(k1+2*k2+2*k3+k4)/6
    evolved = field_from_packed(system.geometry, magnetic+dt*(b1+2*b2+2*b3+b4)/6, field.axial_bias)
    system.fluid.fields(updated)
    # Validate actual sampled density as well as positive point values.
    sampled = system.mapping.scalar(updated[..., 0])
    if np.any(sampled <= 0) or not np.isfinite(sampled).all():
        raise ValueError('Nonpositive evolved reconstructed density; no input modified')
    return updated, evolved, dict(max_residual=max(s.residual for s in stages),
                                  max_iterations=max(s.diagnostics['iterations'] for s in stages))
