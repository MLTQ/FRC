"""Enrich mechanical resolution using existing, well-separated density faces."""
import numpy as np
from .precise_geometry import decimal, geometry_precision
from .precise_search import searchsorted


@geometry_precision
def refine_velocity_faces(density_faces, velocity_faces, rounds, maximum_cells):
    """Keep every old knot and add available middle-half midpoint candidates."""
    if not isinstance(rounds, int) or rounds < 0:
        raise ValueError('Nonnegative integer velocity refinement rounds required')
    if not isinstance(maximum_cells, int) or maximum_cells < 1:
        raise ValueError('Positive declared velocity refinement capacity required')
    knots = np.asarray(velocity_faces, object).copy()
    initial_cells, completed, capacity_limited = len(knots)-1, 0, False
    for _ in range(rounds):
        left, right = knots[:-1], knots[1:]
        targets = (left+right)/decimal(2)
        upper = np.clip(searchsorted(density_faces, targets), 1, len(density_faces)-1)
        lower = upper-1
        # Ties select the lower existing face, matching the documented stage study.
        indices = np.where(abs(density_faces[lower]-targets) <= abs(density_faces[upper]-targets), lower, upper)
        candidate = density_faces[indices]
        eligible = ((candidate >= left+(right-left)/decimal(4)) &
                    (candidate <= right-(right-left)/decimal(4)))
        refined = np.unique(np.r_[knots, candidate[eligible]])
        if len(refined)-1 > maximum_cells:
            capacity_limited = True
            break
        knots = refined
        completed += 1
    return knots, dict(velocity_refinement_requested_rounds=rounds,
        velocity_refinement_completed_rounds=completed,
        velocity_refinement_initial_cells=initial_cells,
        velocity_refinement_added_nodes=len(knots)-1-initial_cells,
        velocity_refinement_capacity=maximum_cells,
        velocity_refinement_capacity_limited=capacity_limited)
