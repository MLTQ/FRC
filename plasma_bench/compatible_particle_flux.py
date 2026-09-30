"""One conservative continuity law and its physical, divergence-commuting flux."""
from dataclasses import dataclass
import numpy as np
from .particle_flux_basis import primitive_basis,difference
from .normal_momentum_faces import radial_faces,axial_faces


@dataclass(frozen=True)
class ParticleFaces:
    radial: np.ndarray
    axial: np.ndarray


class CompatibleParticleFlux:
    """Normal fluxes are qx=r*n*vr and qz=n*vz, with transverse cell means."""
    def __init__(self,density_map):
        self.density_map=density_map;self.geometry=density_map.geometry;self.mesh=density_map.mesh
        geo=self.geometry;mesh=self.mesh
        if geo.order<2*geo.degree+1:
            raise ValueError('Bulk continuity/entropy work requires radial quadrature >= 2*degree+1')
        self.dx=difference(.5*mesh.rf**2);self.dz=difference(mesh.zf,True)
        _,self.hx=primitive_basis(.5*mesh.rf**2,.5*geo.radial_points**2,geo.degree)
        _,self.hz=primitive_basis(mesh.zf,geo.axial_points,geo.degree,periodic=True)
        self.radial_face=radial_faces(mesh,geo.degree,power=1)
        self.axial_face=axial_faces(mesh,geo.degree)

    def validate(self,faces):
        radial=np.asarray(faces.radial,float);axial=np.asarray(faces.axial,float)
        if (radial.shape!=(self.mesh.nr+1,self.mesh.nz) or axial.shape!=(self.mesh.nr,self.mesh.nz)
            or not np.isfinite(radial).all() or not np.isfinite(axial).all()):
            raise ValueError('Finite radial and distinct periodic axial particle faces required')
        if np.any(radial[[0,-1]]!=0):raise ValueError('Zero axis/wall particle flux required')
        return radial,axial

    def from_momentum(self,particle_momentum):
        """Use actual cell means of (n vr,n vz), never analytic point velocities."""
        p=np.asarray(particle_momentum,float)
        if p.shape!=(self.mesh.nr,self.mesh.nz,2) or not np.isfinite(p).all():
            raise ValueError('Finite conservative particle-momentum cell means required')
        radial=self.mesh.rf[:,None]*(self.radial_face@p[...,0]);radial[[0,-1]]=0.
        axial=np.asarray(p[...,1]@self.axial_face.T)
        return ParticleFaces(radial,axial)

    def density_rate(self,faces):
        radial,axial=self.validate(faces)
        return -np.asarray(self.dx@radial+axial@self.dz.T)

    def sample(self,faces):
        radial,axial=self.validate(faces);geo=self.geometry;mapping=self.density_map
        return np.stack((geo.tensor(self.hx,radial,mapping.axial)/geo.r,
                         np.zeros_like(geo.r),geo.tensor(mapping.radial,axial,self.hz)),axis=-1)

    def density_adjoint(self,load):
        """Face coefficients of load dot n_t, with load already volume-weighted."""
        load=self.density_map.validate(load)
        radial=-np.asarray(self.dx.T@load);radial[[0,-1]]=0.
        return ParticleFaces(radial,-np.asarray(load@self.dz))

    def vector_adjoint(self,values):
        """Face coefficients of integral sampled_flux dot values dV."""
        geo=self.geometry;mapping=self.density_map;values=np.asarray(values,float)
        if values.shape!=geo.r.shape+(3,) or not np.isfinite(values).all():
            raise ValueError('Finite physical vector at common quadrature required')
        radial=np.asarray((self.hx.T@geo._unlayout(geo.weight*values[...,0]/geo.r))@mapping.axial)
        axial=np.asarray((mapping.radial.T@geo._unlayout(geo.weight*values[...,2]))@self.hz)
        radial[[0,-1]]=0.
        return ParticleFaces(radial,axial)

    def integrated(self,faces):
        """Actual particle/s through each oriented full annular face."""
        radial,axial=self.validate(faces);mesh=self.mesh
        return 2*np.pi*mesh.dz*radial,2*np.pi*mesh.radial_measure[:,None]*axial
