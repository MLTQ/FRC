"""Shared face fluxes, corner electric fields and cylindrical metric source."""
from dataclasses import dataclass
import numpy as np


@dataclass
class FluxUpdate:
    radial: np.ndarray
    axial: np.ndarray
    electric: np.ndarray
    radial_source: np.ndarray

    def rhs(self, mesh):
        derivative = -np.diff(mesh.rf[:, None, None]*self.radial, axis=0)/mesh.radial_measure[:, None, None]
        derivative -= np.diff(self.axial, axis=1)/mesh.dz
        derivative[..., 1] += self.radial_source
        outward = 2*np.pi*(mesh.radius*mesh.dz*np.sum(self.radial[-1, :, 3])+
                          np.sum(mesh.radial_measure*(self.axial[:, -1, 3]-self.axial[:, 0, 3])))
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
                            np.where(troubled, low.radial_source, self.radial_source))
        # A shared corner affects all four surrounding cells, including diagonals.
        affected = corner[:-1, :-1] | corner[1:, :-1] | corner[:-1, 1:] | corner[1:, 1:]
        return result, affected
