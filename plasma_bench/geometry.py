"""Validated cylindrical geometries and uniform nodal meshes."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Cylinder:
    radius: float = 0.05
    half_length: float = 0.15

    def __post_init__(self):
        if not all(np.isfinite(x) and x > 0 for x in (self.radius, self.half_length)):
            raise ValueError("Cylinder dimensions must be finite and positive")

    def contains(self, r, z):
        tol = 16*np.finfo(float).eps*max(self.radius, self.half_length)
        return (r <= self.radius+tol) & (np.abs(z) <= self.half_length+tol)


@dataclass(frozen=True)
class Grid:
    radius: float = 0.15
    half_length: float = 0.30
    nr: int = 61
    nz: int = 121

    def __post_init__(self):
        Cylinder(self.radius, self.half_length)
        if any(not isinstance(n, int) or isinstance(n, bool) or n < 5 for n in (self.nr, self.nz)):
            raise ValueError("Grid needs at least five integer nodes in each direction")

    @property
    def r(self):
        return np.linspace(0, self.radius, self.nr)

    @property
    def z(self):
        return np.linspace(-self.half_length, self.half_length, self.nz)

    @property
    def dr(self):
        return self.radius / (self.nr - 1)

    @property
    def dz(self):
        return 2 * self.half_length / (self.nz - 1)

    @property
    def mesh(self):
        return np.meshgrid(self.r, self.z, indexing="ij")

    @property
    def shape(self):
        return self.nr, self.nz
