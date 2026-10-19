"""Local reconstruction bounds; never modify conserved cell averages."""
import numpy as np


TRACE_MARGIN = 1e-10  # Relative to each positive cell value, not an SI floor.


def positive_scale(center, proposed_minimum, bound=None):
    """Scale a scalar polynomial about its mean to bound sampled values.

    `bound` (per cell, <= center) replaces the default relative trace margin,
    e.g. a local minimum of neighbouring cell averages.
    """
    bound = TRACE_MARGIN*center if bound is None else np.minimum(bound, (1-1e-12)*center)
    needs = proposed_minimum < bound
    factor = np.ones_like(center)
    np.divide(center-bound, center-proposed_minimum, out=factor, where=needs)
    return np.where(needs, np.clip(factor, 0, 1)*(1-32*np.finfo(float).eps), 1.)


def positive_blend(positive, proposed, center, bound=None):
    """Largest weight of proposed in a convex blend with an admissible trace."""
    bound = TRACE_MARGIN*center if bound is None else np.minimum(bound, (1-1e-12)*center)
    needs = proposed < bound
    weight = np.ones_like(center)
    np.divide(positive-bound, positive-proposed, out=weight, where=needs)
    return np.where(needs, np.clip(weight, 0, 1)*(1-32*np.finfo(float).eps), 1.)


def neighbour_minimum(values):
    """Minimum over each cell and its 3x3 neighbours (edge-replicated at array edges)."""
    padded = np.pad(values, ((1, 1), (1, 1)), mode="edge")
    rows, cols = values.shape
    return np.min([padded[a:a+rows, b:b+cols] for a in range(3) for b in range(3)], axis=0)


def limiter_statistics(density_scale, thermal_scale, pressure_weight):
    return {"density_limited_cells": int(np.count_nonzero(density_scale < 1)),
            "thermal_slope_limited_cells": int(np.count_nonzero(thermal_scale < 1)),
            "pressure_blended_cells": int(np.count_nonzero(pressure_weight < 1)),
            "minimum_pressure_blend_weight": float(np.min(pressure_weight))}


def empty_statistics():
    return {"density_limited_cells": 0, "thermal_slope_limited_cells": 0,
            "pressure_blended_cells": 0, "minimum_pressure_blend_weight": 1.}
