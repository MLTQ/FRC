"""Immutable conservative density reconstruction and its physical adjoint."""
import numpy as np
from .compatible_velocity import radial_moments, axial_moments


class CompatibleDensity:
    """Number-density cell averages; rates come from the caller's continuity law."""
    def __init__(self, geometry, cell_means):
        self.geometry=geometry;self.mesh=geometry.mesh
        self.radial=radial_moments(self.mesh,geometry.radial_points,geometry.degree)
        self.axial=axial_moments(self.mesh,geometry.axial_points,geometry.degree)
        self.cell_means=self.validate(cell_means).copy()
        self.samples=self.sample(self.cell_means)
        if np.any(self.cell_means<=0) or np.any(self.samples<=0):
            raise ValueError('Positive cell and reconstructed density required; no floor or fallback')
        self.cell_means.setflags(write=False);self.samples.setflags(write=False)

    def validate(self, values):
        values=np.asarray(values,float)
        if values.shape!=(self.mesh.nr,self.mesh.nz) or not np.isfinite(values).all():
            raise ValueError('Finite number-density cell averages or their rates required')
        return values

    def sample(self, values):
        with np.errstate(over='ignore',invalid='ignore'):
            result=self.geometry.tensor(self.radial,self.validate(values),self.axial)
        if not np.isfinite(result).all():
            raise ValueError('Density reconstruction or rate exceeds finite numerical range')
        return result

    def adjoint(self, values):
        """R.T W values, including physical volume; not divided by cell volumes."""
        geo=self.geometry;values=np.asarray(values,float)
        if values.shape!=geo.r.shape or not np.isfinite(values).all():
            raise ValueError('Finite scalar on the common physical quadrature required')
        weighted=geo._unlayout(geo.weight*values)
        return np.asarray((self.radial.T@weighted)@self.axial)
