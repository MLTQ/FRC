"""Explicit target-relative spacing for a newly constructed velocity mesh."""
import numpy as np
from .precise_geometry import decimal, geometry_precision


@geometry_precision
def select_velocity_indices(faces, targets, minimum_target_fraction=0.):
    if not np.isfinite(minimum_target_fraction) or not 0 <= minimum_target_fraction <= 1:
        raise ValueError('Target spacing fraction must lie in [0,1]')
    nearest = np.argmin(abs(faces[:, None]-targets[None, :]), axis=0)
    candidates = np.unique(nearest)
    if not minimum_target_fraction:
        return candidates
    # Each candidate represents the target it approximates most closely.
    ranks = [min(np.flatnonzero(nearest == c), key=lambda j: abs(faces[c]-targets[j]))
             for c in candidates]
    chosen = list(zip(candidates.tolist(), ranks))
    fraction = decimal(minimum_target_fraction)
    while len(chosen) > 2:
        for k, ((a, i), (b, j)) in enumerate(zip(chosen[:-1], chosen[1:])):
            if faces[b]-faces[a] >= fraction*(targets[j]-targets[i]):
                continue
            # Physical support endpoints are mandatory. No old velocity knot
            # enters this routine: it constructs only a previously absent space.
            if k == 0:
                remove = k+1
            elif k+1 == len(chosen)-1:
                remove = k
            else:
                remove = k if abs(faces[a]-targets[i]) > abs(faces[b]-targets[j]) else k+1
            chosen.pop(remove)
            break
        else:
            break
    return np.asarray([c for c, _ in chosen], int)
