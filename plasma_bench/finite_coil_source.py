"""Finite circular-coil source and affine induced/source magnetic fields."""
from dataclasses import dataclass
import numpy as np
from .electromagnetics import loop_field


class FiniteCoilSource:
    def __init__(self,geometry,coils):
        self.geometry=geometry;self.coils=tuple(coils)
        self.psi_unit,self.b_unit=self.evaluate(geometry.r,geometry.z)
        self.a_unit=np.zeros_like(self.b_unit)
        self.a_unit[...,1]=self.psi_unit/geometry.r

    def evaluate(self,r,z):
        r,z=np.broadcast_arrays(r,z);psi=np.zeros_like(r,dtype=float);b=np.zeros(r.shape+(3,))
        for radius,center,turns in self.coils:
            if radius<=self.geometry.mesh.radius and abs(center)<=self.geometry.mesh.half_length:
                raise ValueError('Prescribed source filament must be outside plasma')
            p,br,bz=loop_field(r,z,radius,center)
            psi+=turns*p;b[...,0]+=turns*br;b[...,2]+=turns*bz
        return psi,b

    def linkage(self,current):
        g=self.geometry
        return float(np.sum(g.weight*np.sum(g.sample_electric(current)*self.a_unit,axis=-1)))


@dataclass(frozen=True)
class FiniteCoilMagnetic:
    source: FiniteCoilSource
    psi_increment: np.ndarray
    g: np.ndarray
    source_amplitude: float

    def packed_increment(self,geometry):
        if geometry is not self.source.geometry or not np.isfinite(self.source_amplitude):
            raise ValueError('Same finite source geometry and finite current required')
        return geometry.pack(self.psi_increment,self.g)

    def samples(self,geometry):
        self.packed_increment(geometry)
        return geometry.magnetic_samples(self.psi_increment,self.g)+self.source_amplitude*self.source.b_unit
