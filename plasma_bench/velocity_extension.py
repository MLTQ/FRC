"""Retain old velocity nodes and explicitly extend their edge polynomials."""
import numpy as np
from .precise_geometry import geometry_precision


@geometry_precision
def inherited_velocity_faces(faces, old, *, extend=False):
    inside = old[(old >= faces[0]) & (old <= faces[-1])]
    if not len(inside):
        return None, False
    if not extend or len(inside) < 2:
        return np.unique(np.r_[faces[0], inside, faces[-1]]), False
    left = faces[0] < inside[0] and inside[0]-faces[0] <= inside[1]-inside[0]
    right = faces[-1] > inside[-1] and faces[-1]-inside[-1] <= inside[-1]-inside[-2]
    selected = np.r_[[] if left else [faces[0]], inside, [] if right else [faces[-1]]]
    return np.unique(selected), bool(left or right)
