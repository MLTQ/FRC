"""Checked implicit midpoint for the paired electron pressure split component."""
from dataclasses import dataclass
import numpy as np
from .thermal_means import logarithmic_mean


class ElectronConvergenceError(RuntimeError):
    """No input state has been committed."""


@dataclass(frozen=True)
class ElectronStep:
    g: np.ndarray
    entropy_density: np.ndarray
    iterations: int
    residual: float


def step(model, g, entropy_density, dt, *, tolerance=1e-11, max_iterations=80):
    if not np.isfinite([dt, tolerance]).all() or dt <= 0 or not 0 < tolerance < 1 or type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('Positive finite time/tolerance and positive iteration budget required')
    oldg, oldu = np.asarray(g, float).copy(), np.asarray(entropy_density, float).copy()
    dg, du, _, _ = model.rhs(oldg, oldu)
    gscale = max(float(abs(oldg).max()), float(dt*abs(dg).max()), np.finfo(float).tiny)
    uscale = max(float(abs(oldu).max()), float(dt*abs(du).max()), np.finfo(float).tiny)
    trialg, trialu = oldg.copy(), oldu.copy()
    for iteration in range(1, max_iterations+1):
        dg, du, _, _ = model.rhs(.5*(oldg+trialg), logarithmic_mean(oldu, trialu))
        nextg, nextu = oldg+dt*dg, oldu+dt*du
        if not np.isfinite(nextg).all() or not np.isfinite(nextu).all() or np.any(nextu <= 0):
            raise ElectronConvergenceError('Nonpositive/nonfinite electron trial; reduce the step, no state repair applied')
        residual = max(float(abs(nextg-trialg).max())/gscale, float(abs(nextu-trialu).max())/uscale)
        trialg, trialu = nextg, nextu
        if residual <= tolerance:
            # Check the actual midpoint equations at the returned state.
            dg, du, _, _ = model.rhs(.5*(oldg+trialg), logarithmic_mean(oldu, trialu))
            residual = max(float(abs(trialg-oldg-dt*dg).max())/gscale, float(abs(trialu-oldu-dt*du).max())/uscale)
            if residual <= tolerance:
                return ElectronStep(trialg, trialu, iteration, residual)
    raise ElectronConvergenceError(f'Electron midpoint failed after {max_iterations} iterations; residual {residual:.6g}')
