"""Conservative overlap of primal mass and its flux onto nodal dual volumes."""
import numpy as np
from .dual_current import DualCurrent


class DualMass:
    def __init__(self, mesh):
        self.mesh = mesh
        self.geometry = DualCurrent(mesh)
        self.volumes = self.geometry.volumes
        edges = self.geometry.dual_x_faces
        self.centers_x = .5*(edges[:-1]+edges[1:])[:, None]

    def density(self, primal):
        """Overlap averages; each primal cell contributes four quarter masses."""
        value = np.asarray(primal, float)
        if value.shape != (self.mesh.nr, self.mesh.nz) or not np.isfinite(value).all():
            raise ValueError('Finite primal cell density or density rate required')
        quarter = .25*self.mesh.volumes*value
        pair = quarter+np.roll(quarter, 1, axis=1)
        mass = np.zeros_like(self.volumes)
        mass[:-1] += pair
        mass[1:] += pair
        return mass/self.volumes

    def primal_average(self, dual):
        """Cell average of the nodal Q1 scalar; adjoint to the overlap map."""
        value = self.geometry.scalar(dual)
        pair = value[:-1]+value[1:]
        return .25*(pair+np.roll(pair, -1, axis=1))

    def fluxes(self, radial, axial):
        """Integrate the primal RT0 flux over each dual face, in the same units."""
        radial, axial = np.asarray(radial, float), np.asarray(axial, float)
        mesh = self.mesh
        if radial.shape != (mesh.nr+1, mesh.nz) or axial.shape != (mesh.nr, mesh.nz+1):
            raise ValueError('Integrated primal radial and duplicate-periodic axial face fluxes required')
        if not np.isfinite(radial).all() or not np.isfinite(axial).all() or not np.array_equal(axial[:, 0], axial[:, -1]):
            raise ValueError('Finite fluxes with identical periodic endpoint values required')
        center = .5*(radial[:-1]+radial[1:])
        dual_r = np.empty((mesh.nr+2, mesh.nz))
        dual_r[1:-1] = .5*(center+np.roll(center, 1, axis=1))
        dual_r[[0, -1]] = .5*(radial[[0, -1]]+np.roll(radial[[0, -1]], 1, axis=1))
        center = .25*(axial[:, :-1]+axial[:, 1:])
        dual_z = np.zeros_like(self.volumes)
        dual_z[:-1] += center
        dual_z[1:] += center
        return dual_r, dual_z
