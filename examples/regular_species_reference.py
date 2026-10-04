"""Analytic original species pressure/entropy rates, independent of the pair."""
import numpy as np
from scipy.constants import elementary_charge as E
from .regular_coupling_reference import initialize as initialize_base


def entropy(ref, r, z, species):
    y = (r/ref.radius)**2
    phase = np.pi*z/ref.half_length
    return (np.log(2. if species == 0 else 1.5)
            + (.13 if species == 0 else -.09)*y*np.sin(phase)
            + .07*np.cos(2*phase))


def initialize(spec):
    ref, fluid, state, model, field = initialize_base(spec)
    for k in range(2):
        state[..., k+4] = state[..., 0]*entropy(ref, fluid.r, fluid.z, k)
    return ref, fluid, state, model, field


def exact(ref, r, z, gammas, n0, fractions):
    """Original grad(p)/(e n) and -a J dot grad(s)/e; analytic derivatives."""
    r, z = np.broadcast_arrays(r, z)
    f = ref.fields(r, z)
    y = (r/ref.radius)**2
    wave = np.pi/ref.half_length
    phase = wave*z
    pressure, rates = [], []
    for k, (gamma, a) in enumerate(zip(gammas, fractions)):
        s = entropy(ref, r, z, k)
        b = .13 if k == 0 else -.09
        ds = np.stack((b*2*r/ref.radius**2*np.sin(phase), np.zeros_like(r),
                       b*y*wave*np.cos(phase)-.14*wave*np.sin(2*phase)), axis=-1)
        p = E*np.exp(s)*f['n']*(f['n']/n0)**(gamma-1)
        dp = p[..., None]*(gamma*f['gradn']/f['n'][..., None]+ds)
        pressure.append(a*dp/(E*f['n'][..., None]))
        rates.append(-a/E*np.sum(f['current']*ds, axis=-1))
    return dict(pressure=np.stack(pressure), entropy_rate=np.stack(rates))
