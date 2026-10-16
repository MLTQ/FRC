"""Conservative linear Euler traces with a common admissibility slope factor."""
import numpy as np


def thermal(conserved, mass_ratio):
    number = conserved[..., 0]
    momentum = conserved[..., 1:4]
    velocity = np.divide(momentum, number[..., None], out=np.zeros_like(momentum), where=number[..., None] > 0)/mass_ratio[:, None]
    return conserved[..., 4]-.5*mass_ratio*number*np.sum(velocity**2, axis=-1), velocity


def reconstruct(conserved, faces, mass_ratio, gamma=5/3):
    """Cell-mean-conservative endpoint traces, without modifying cell states."""
    width = np.diff(faces); centre = .5*(faces[:-1]+faces[1:])
    gradient = np.diff(conserved, axis=0)/np.diff(centre)[:, None, None]
    slope = np.zeros_like(conserved); slope[0] = gradient[0]; slope[-1] = gradient[-1]
    slope[1:-1] = np.where(gradient[:-1]*gradient[1:] > 0,
        np.sign(gradient[:-1])*np.minimum(abs(gradient[:-1]), abs(gradient[1:])), 0.)
    delta = .5*width[:, None, None]*slope

    def admissible(q):
        internal, _ = thermal(q, mass_ratio)
        empty = (q[..., 0] == 0)&np.all(q[..., 1:] == 0, axis=-1)
        return empty|((q[..., 0] > 0)&(internal >= 0)&np.isfinite(q).all(axis=-1))

    high = np.ones(conserved.shape[:2]); low = np.zeros_like(high)
    unlimited = admissible(conserved-delta)&admissible(conserved+delta)
    low[unlimited] = 1.
    # The Euler admissible set is convex. Scale the entire conservative slope,
    # retaining the cell mean and cross-component kinetic-energy consistency.
    if not np.all(unlimited):
        for _ in range(48):
            middle = .5*(low+high); d = middle[..., None]*delta
            good = admissible(conserved-d)&admissible(conserved+d)
            low = np.where(good, middle, low); high = np.where(good, high, middle)
    delta *= low[..., None]

    def primitive(q):
        internal, velocity = thermal(q, mass_ratio)
        return np.concatenate((q[..., :1], ((gamma-1)*internal)[..., None], velocity), axis=-1)

    return primitive(conserved-delta), primitive(conserved+delta), low
