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
    toroidal_axis: np.ndarray | None = None   # 2.5D: E_z on the axis (nz,), the r=0 term of r F_r = -E_z for B_theta/r
    toroidal_boundary: tuple | None = None    # 2.5D, linked vacuum r B_theta only: boundary B_theta Poynting flux
                                              # (wall (nz,), bottom (nr,), top (nr,), inner (nz,) or None)

    def rhs(self, mesh):
        radial, axial = self.radial, self.axial
        if self.boundary_poynting is not None:
            radial, axial = radial.copy(), axial.copy()
            wall, bottom, top, inner = self.boundary_poynting.energy_fluxes(mesh, self.electric)
            if self.toroidal_boundary is not None:   # the poloidal variational flux plus the vacuum B_theta part
                tw, tb, tt, ti = self.toroidal_boundary
                wall, bottom, top = wall+tw, bottom+tb, top+tt
                inner = None if inner is None else inner+ti
            radial[-1, :, 3] = wall
            axial[:, 0, 3], axial[:, -1, 3] = bottom, top
            if inner is not None:
                radial[0, :, 3] = inner
        derivative = -np.diff(mesh.rf[:, None, None]*radial, axis=0)/mesh.radial_measure[:, None, None]
        derivative -= np.diff(axial, axis=1)/mesh.dz
        derivative[..., 1] += self.radial_source
        if self.toroidal_axis is not None:   # dg/dt = (E_z(r_1) - E_z(0))/integral(r dr) in the axis cells
            derivative[0, :, 5] -= self.toroidal_axis/mesh.radial_measure[0]
        outward = 2*np.pi*(mesh.radius*mesh.dz*np.sum(radial[-1, :, 3])+
                          np.sum(mesh.radial_measure*(axial[:, -1, 3]-axial[:, 0, 3])))
        if mesh.inner_radius > 0:   # energy leaving through an inner wall (toward the axis)
            outward -= 2*np.pi*mesh.inner_radius*mesh.dz*np.sum(radial[0, :, 3])
        if self.toroidal_boundary is not None:   # 4th entry: B_theta Poynting power INTO the gas (part of -outward)
            tw, tb, tt, ti = self.toroidal_boundary
            inflow = -2*np.pi*(mesh.radius*mesh.dz*np.sum(tw)+np.sum(mesh.radial_measure*(tt-tb)))
            if ti is not None:
                inflow += 2*np.pi*mesh.inner_radius*mesh.dz*np.sum(ti)
            return derivative, -mesh.rf[:, None]*self.electric, float(outward), float(inflow)
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
                            np.where(troubled, low.radial_source, self.radial_source), self.boundary_poynting,
                            None if self.toroidal_axis is None else np.where(troubled[0], low.toroidal_axis, self.toroidal_axis),
                            None if self.toroidal_boundary is None else tuple(
                                None if mine is None else np.where(mask, base, mine) for mine, base, mask in
                                zip(self.toroidal_boundary, low.toroidal_boundary,
                                    (troubled[-1], troubled[:, 0], troubled[:, -1], troubled[0]))))
        # A shared corner affects all four surrounding cells, including diagonals.
        affected = corner[:-1, :-1] | corner[1:, :-1] | corner[:-1, 1:] | corner[1:, 1:]
        return result, affected
