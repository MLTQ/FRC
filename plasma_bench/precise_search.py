"""Exact precise-coordinate search within monotone float brackets."""
import numpy as np


def searchsorted(edges, points, side='left'):
    if side not in ('left', 'right'):
        raise ValueError('Search side must be left or right')
    query = np.asarray(points, object)
    if len(edges) < 512 or query.size < 128:
        return np.searchsorted(edges, query, side=side)
    flat = query.ravel()
    high_edges, high_points = np.asarray(edges, float), np.asarray(flat, float)
    # Strictly different rounded values determine exact ordering. Only equal
    # rounded keys need exact comparisons; all ties retain their precise values.
    lo = np.searchsorted(high_edges, high_points, side='left')
    hi = np.searchsorted(high_edges, high_points, side='right')
    active = np.flatnonzero(lo < hi)
    while active.size:
        mid = lo[active]+(hi[active]-lo[active])//2
        below = edges[mid] < flat[active] if side == 'left' else edges[mid] <= flat[active]
        lower, upper = active[below], active[~below]
        lo[lower] = mid[below]+1
        hi[upper] = mid[~below]
        active = active[lo[active] < hi[active]]
    return lo.reshape(query.shape)

