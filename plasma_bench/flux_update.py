"""Shared face fluxes, corner electric fields and cylindrical metric source."""
from dataclasses import dataclass
import numpy as np


@dataclass
class FluxUpdate:
    radial: np.ndarray
    axial: np.ndarray
    electric: np.ndarray
    radial_source: np.ndarray
    boundary_poynting: object | None = None

    def rhs(self, mesh):
        radial, axial = self.radial, self.axial
        if self.boundary_poynting is not None:
            radial, axial = radial.copy(), axial.copy()
            wall, bottom, top = self.boundary_poynting.energy_fluxes(mesh, self.electric)
            radial[-1, :, 3] = wall
            axial[:, 0, 3], axial[:, -1, 3] = bottom, top
        derivative = -np.diff(mesh.rf[:, None, None]*radial, axis=0)/mesh.radial_measure[:, None, None]
        derivative -= np.diff(axial, axis=1)/mesh.dz
        derivative[..., 1] += self.radial_source
        outward = 2*np.pi*(mesh.radius*mesh.dz*np.sum(radial[-1, :, 3])+
                          np.sum(mesh.radial_measure*(axial[:, -1, 3]-axial[:, 0, 3])))
        return derivative, -mesh.rf[:, None]*self.electric, float(outward)

    def local_fallback(self, low, troubled):
        """Replace every incident flux/source of troubled cells by its baseline."""
        radial = np.zeros(self.radial.shape[:2], dtype=bool)
        axial = np.zeros(self.axial.shape[:2], dtype=bool)
        corner = np.zeros(self.electric.shape, dtype=bool)
        radial[:-1] |= troubled; radial[1:] |= troubled
        axial[:, :-1] |= troubled; axial[:, 1:] |= troubled
        for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
            corner[i:i+troubled.shape[0], j:j+troubled.shape[1]] |= troubled
        result = FluxUpdate(np.where(radial[..., None], low.radial, self.radial),
                            np.where(axial[..., None], low.axial, self.axial),
                            np.where(corner, low.electric, self.electric),
                            np.where(troubled, low.radial_source, self.radial_source), self.boundary_poynting)
        # A shared corner affects all four surrounding cells, including diagonals.
        affected = corner[:-1, :-1] | corner[1:, :-1] | corner[:-1, 1:] | corner[1:, 1:]
        return result, affected
