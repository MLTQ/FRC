"""Annular finite volumes and a divergence-free corner-flux representation."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MHDMesh:
    radius: float
    half_length: float
    nr: int
    nz: int

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in (self.radius, self.half_length)):
            raise ValueError("Mesh dimensions must be finite and positive")
        if any(not isinstance(n, int) or isinstance(n, bool) or n < 4 for n in (self.nr, self.nz)):
            raise ValueError("MHD mesh needs at least four cells per direction")

    @property
    def dr(self):
        return self.radius/self.nr

    @property
    def dz(self):
        return 2*self.half_length/self.nz

    @property
    def rf(self):
        return np.linspace(0, self.radius, self.nr+1)

    @property
    def zf(self):
        return np.linspace(-self.half_length, self.half_length, self.nz+1)

    @property
    def r(self):
        return (self.rf[:-1]+self.rf[1:])/2

    @property
    def z(self):
        return (self.zf[:-1]+self.zf[1:])/2

    @property
    def centers(self):
        return np.meshgrid(self.r, self.z, indexing="ij")

    @property
    def corners(self):
        return np.meshgrid(self.rf, self.zf, indexing="ij")

    @property
    def radial_measure(self):
        return np.diff(self.rf**2)/2

    @property
    def radial_centroid(self):
        return np.diff(self.rf**3)/(3*self.radial_measure)

    @property
    def volumes(self):
        return np.broadcast_to(2*np.pi*self.radial_measure[:, None]*self.dz, (self.nr, self.nz))

    def face_fields(self, psi):
        psi = np.asarray(psi)
        if psi.shape != (self.nr+1, self.nz+1):
            raise ValueError("psi must be defined on cell corners")
        denominator = self.rf[:, None]*self.dz
        br = np.divide(-np.diff(psi, axis=1), denominator,
                       out=np.zeros((self.nr+1, self.nz)), where=denominator > 0)
        bz = np.diff(psi, axis=0)/self.radial_measure[:, None]
        return br, bz

    def cell_fields(self, psi):
        br, bz = self.face_fields(psi)
        weight = ((self.radial_centroid-self.rf[:-1])/self.dr)[:, None]
        return (1-weight)*br[:-1]+weight*br[1:], (bz[:, :-1]+bz[:, 1:])/2

    def divergence(self, psi):
        br, bz = self.face_fields(psi)
        radial = np.diff(self.rf[:, None]*br, axis=0)/self.radial_measure[:, None]
        return radial+np.diff(bz, axis=1)/self.dz
