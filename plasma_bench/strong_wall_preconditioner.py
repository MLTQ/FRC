"""Restrict a fixed-background approximate inverse to the physical wall space."""
import numpy as np
from .plasma_preconditioner import PlasmaPreconditioner
from .regular_current import RegularCurrent


def magnetic_metric(coordinates, state):
    """Full-coordinate physical magnetic plus current-inertia Hessian."""
    c = coordinates; geo = c.geometry
    model = RegularCurrent(c.system.mapping, state[..., 0], c.system.ion_mass)
    currents = []
    for k, index in enumerate(c.magnetic_indices):
        b = np.zeros(geo.hessian.shape[0]); b[index] = c.magnetic_scale[k]
        currents.append(model.current(b))
    currents = np.column_stack(currents)
    metric = geo.hessian[np.ix_(c.magnetic_indices, c.magnetic_indices)].toarray()
    metric *= c.magnetic_scale[:, None]*c.magnetic_scale[None, :]
    return metric+currents.T@np.column_stack([model.inertia@j for j in currents.T])


def magnetic_left(coordinates, state):
    q = coordinates.magnetic_basis
    metric = magnetic_metric(coordinates.full, state)
    return np.linalg.solve(q.T@metric@q, q.T@metric)


class StrongWallPreconditioner(PlasmaPreconditioner):
    def __init__(self, coordinates):
        self.coordinates = c = coordinates; original = PlasmaPreconditioner(c.full)
        keep = c.fluid_keep; q = c.magnetic_basis
        left = magnetic_left(c, c.background())
        self.fluid = original.fluid[keep][:, keep].tocsc()
        self.to_fluid = original.to_fluid[keep]@q
        self.to_magnetic = left@original.to_magnetic[:, keep]
        self.magnetic = left@original.magnetic@q
        self._dt = None
