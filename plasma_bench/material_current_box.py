"""Bounded pressure trace for the actual implicit current substep."""
import numpy as np
from scipy.optimize import minimize


def box_minimum(matrix, linear, limits, guess):
    """Solve the convex current quadratic, retaining compatible null gradients."""
    def objective(y):
        product = matrix@y
        return float(.5*y@product-linear@y), product-linear
    fit = minimize(objective, np.clip(guess, -limits, limits), jac=True,
        bounds=list(zip(-limits, limits)), method='L-BFGS-B',
        options={'gtol': 2e-12, 'ftol': 0., 'maxiter': 1000, 'maxls': 40})
    y = fit.x.copy()
    working = np.zeros(len(y), dtype=int)
    working[y <= -limits+1e-14] = -1
    working[y >= limits-1e-14] = 1
    y[working != 0] = working[working != 0]*limits[working != 0]
    for iteration in range(4*len(y)+20):
        gradient = matrix@y-linear
        error = float(abs(y-np.clip(y-gradient, -limits, limits)).max())
        if error < 2e-12:
            return y
        free = np.flatnonzero(working == 0)
        if len(free):
            block = matrix[np.ix_(free, free)]; scale = np.sqrt(np.diag(block))
            # The curl-compatible metric is positive semidefinite. Preserve its
            # null gradients; adding a diagonal penalty would change the model.
            delta = np.linalg.lstsq(block/scale[:, None]/scale[None, :],
                -gradient[free]/scale, rcond=1e-13)[0]/scale
            alpha = 1.; hit = None; side = 0
            for direction in (-1, 1):
                moving = np.flatnonzero(direction*delta > 0)
                if len(moving):
                    ratios = (direction*limits[free[moving]]-y[free[moving]])/delta[moving]
                    k = int(ratios.argmin())
                    if ratios[k] < alpha:
                        alpha = float(ratios[k]); hit = int(free[moving[k]]); side = direction
            y[free] += max(0., alpha)*delta
            if hit is not None:
                y[hit] = side*limits[hit]; working[hit] = side
                continue
        gradient = matrix@y-linear
        violation = working*gradient
        if np.max(violation) > 2e-12:
            # Release one wrong-sign constraint only after the current free
            # equations settle. Releasing every bound on every iteration can
            # create a zero-length outward Newton step at a newly hit bound.
            working[int(violation.argmax())] = 0
        elif np.max(abs(y-np.clip(y-gradient, -limits, limits))) < 2e-12:
            return y
        else:
            break
    raise RuntimeError(f'Upwind midpoint box equations did not converge: {error:g}, {fit.message}')
