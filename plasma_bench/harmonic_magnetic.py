"""Lossless affine storage of a uniform axial harmonic and magnetic increments."""
import numpy as np
from scipy.constants import mu_0 as MU0


class HarmonicMagnetic:
    """A bounded field representation, not a fitted equilibrium or time solver."""
    def __init__(self,mesh,psi_increment,g,axial_bias):
        psi=np.asarray(psi_increment,float);g=np.asarray(g,float)
        if psi.shape!=(mesh.nr+1,mesh.nz+1) or g.shape!=(mesh.nr,mesh.nz):raise ValueError('Complete magnetic increment shapes required')
        if not np.isfinite(psi).all() or not np.isfinite(g).all() or not np.isfinite(axial_bias):raise ValueError('Finite magnetic field required')
        if np.any(psi[[0,-1]]!=0) or np.any(psi[:,0]!=psi[:,-1]):raise ValueError('Increment must vanish at fixed axis/wall and be periodic')
        self.mesh=mesh;self.psi_increment=psi.copy();self.g=g.copy();self.axial_bias=float(axial_bias)
        self.psi_increment.setflags(write=False);self.g.setflags(write=False)

    def packed_increment(self,geometry):
        if geometry.mesh!=self.mesh:raise ValueError('Matching field geometry required')
        return geometry.pack(self.psi_increment,self.g)

    def current(self,model):return model.current(self.packed_increment(model.geometry))

    def samples(self,geometry):
        self.packed_increment(geometry)
        value=geometry.magnetic_samples(self.psi_increment,self.g)
        value[...,2]+=self.axial_bias
        return value

    def energy_components(self,model):
        value=self.packed_increment(model.geometry)
        return dict(harmonic=self.axial_bias**2*float(self.mesh.volumes.sum())/(2*MU0),
                    increment=float(.5*value@(model.hessian@value)),global_cross=0.)

    def energy_density(self,geometry):
        value=self.samples(geometry)
        return np.sum(geometry.weight*np.sum(value*value,axis=-1),axis=(2,3))/(2*MU0*self.mesh.volumes)

    def rounded_total(self,geometry):
        """Explicitly lossy interchange form; never used to compute current."""
        psi=self.psi_increment+.5*self.axial_bias*self.mesh.corners[0]**2
        return geometry.pack(psi,self.g)
