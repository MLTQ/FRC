"""Retain material interfaces in the mechanical space at explicit density contrast."""
import numpy as np
from .precise_geometry import geometry_precision, widths


@geometry_precision
def refine_density_jumps(density_faces, relative_numbers, velocity_faces, minimum_contrast):
    """Add existing density-jump faces; never remove material or inherited knots."""
    faces, knots = np.asarray(density_faces, object), np.asarray(velocity_faces, object)
    number = np.asarray(relative_numbers, float)
    contrast = float(minimum_contrast)
    if (faces.ndim != 1 or number.shape != (len(faces)-1,) or
            not np.isfinite(number).all() or np.any(number <= 0) or
            not np.isfinite(contrast) or contrast <= 1):
        raise ValueError('Positive finite cell amounts and explicit density contrast>1 required')
    h = widths(faces)
    if (np.any(h <= 0) or knots.ndim != 1 or len(knots) < 2 or
            np.any(widths(knots) <= 0) or not np.isin(knots, faces).all()):
        raise ValueError('Ordered density faces and existing nested velocity knots required')
    log_density = np.log(number)-np.log(h)
    log_jumps = abs(np.diff(log_density))
    selected = faces[1:-1][log_jumps >= np.log(contrast)]
    result = np.unique(np.r_[knots, selected])
    return result, dict(density_jump_minimum_contrast=contrast,
        density_jump_selected_faces=len(selected), velocity_refinement_initial_cells=len(knots)-1,
        velocity_refinement_added_nodes=len(result)-len(knots),
        maximum_adjacent_log10_density_jump=float(log_jumps.max(initial=0.)/np.log(10.)))
