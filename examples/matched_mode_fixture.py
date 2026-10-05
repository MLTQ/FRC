"""Independent contact mode on explicitly separate fluid and field meshes."""
import numpy as np
from .strong_wall_reference import initialize
from .continuum_contact_mode import ContinuumContactMode
from .finite_mass_manufactured import cell_average


def fixture(spec):
    """Return continuum reference, unperturbed state and analytic perturbation."""
    mode = ContinuumContactMode()
    if spec.get('kind', 'contact') != 'contact': raise ValueError('Stationary contact mode required')
    _, d, u, system, b = initialize(dict(spec, kind='contact')); u[..., 1:4] = 0.
    target = mode.fields(d.r, d.z)['physical']; du = np.zeros_like(u)
    du[..., 0] = target[..., 0]
    np.divide(u[..., :1]*target[..., 1:3], d.r[..., None], out=du[..., 1:3], where=d.r[..., None] > 0)
    du[-1, -1, :, 1] = 0.; du[..., 3] = u[..., 0]*target[..., 3]
    du[..., 4:] = u[..., :1]*target[..., 4:6]+u[..., 4:]/u[..., :1]*du[..., :1]
    mesh = system.geometry.mesh
    psi = mode.fields(*mesh.corners)['psi']; psi[[0, -1]] = 0.; psi[:, -1] = psi[:, 0]
    g = cell_average(mesh, lambda r, z: mode.fields(r, z)['physical'][..., 7]/r, order=16)
    g = system.geometry.radial_projection@g
    return mode, d, u, system, b, du, system.geometry.pack(psi, g)
