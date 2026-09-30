"""Physical kinetic energy and velocity from true particle momentum means."""
import numpy as np
from .compatible_velocity import radial_moments
from .momentum_derivatives import radial_derivative,axial_derivative


class ConservativeKinematics:
    """Momentum order: (n vr, n vz, n r vtheta); physical vectors: (r,theta,z)."""
    def __init__(self,geometry,total_mass):
        if not np.isfinite(total_mass) or total_mass<=0:raise ValueError('Positive total species mass required')
        self.geometry=geometry;self.mesh=geometry.mesh;self.mass=total_mass
        self.radial={p:radial_moments(self.mesh,geometry.radial_points,geometry.degree,power=p) for p in (1,2)}
        self._derivatives=None

    def validate(self,density,momentum):
        if density.geometry is not self.geometry:raise ValueError('One common physical density/kinetic geometry required')
        p=np.asarray(momentum,float)
        if p.shape!=(self.mesh.nr,self.mesh.nz,3) or not np.isfinite(p).all():
            raise ValueError('Finite conservative particle momentum and angular momentum required')
        return p

    def sampled_momentum(self,density,momentum):
        p=self.validate(density,momentum);g=self.geometry
        return np.stack((g.tensor(self.radial[1],p[...,0],density.axial),
            g.tensor(self.radial[2],p[...,2],density.axial)/g.r,
            g.tensor(density.radial,p[...,1],density.axial)),axis=-1)

    def fields(self,density,momentum):
        p=self.sampled_momentum(density,momentum);v=p/density.samples[...,None]
        energy=.5*self.mass*np.sum(p*v,axis=-1)
        if not np.isfinite(v).all() or not np.isfinite(energy).all():raise ValueError('Kinetic state exceeds finite range')
        return dict(momentum=p,velocity=v,kinetic=energy)

    def average(self,values):
        return np.sum(self.geometry.weight*values,axis=(2,3))/self.mesh.volumes

    def momentum_adjoint(self,density,values):
        g=self.geometry;v=np.asarray(values,float)
        if v.shape!=g.r.shape+(3,) or not np.isfinite(v).all():raise ValueError('Finite physical vector required')
        def pull(radial,value):return np.asarray(radial.T@g._unlayout(g.weight*value)@density.axial)
        return np.stack((pull(self.radial[1],v[...,0]),pull(density.radial,v[...,2]),
                         pull(self.radial[2],v[...,1]/g.r)),axis=-1)

    def energy_gradients(self,density,momentum):
        v=self.fields(density,momentum)['velocity']
        return dict(momentum=self.mass*self.momentum_adjoint(density,v),
            density=density.adjoint(-.5*self.mass*np.sum(v*v,axis=-1)))

    def energy_rate(self,density,momentum,density_rate,momentum_rate):
        v=self.fields(density,momentum)['velocity'];pt=self.sampled_momentum(density,momentum_rate)
        nt=density.sample(density_rate)
        changing_momentum=self.mass*np.sum(v*pt,axis=-1)
        changing_density=-.5*self.mass*np.sum(v*v,axis=-1)*nt
        return dict(momentum=self.average(changing_momentum),density=self.average(changing_density),
                    total=self.average(changing_momentum+changing_density))

    def gradients(self,density,momentum):
        p=self.validate(density,momentum);g=self.geometry;n=density.samples
        if self._derivatives is None:
            self._derivatives=({s:radial_derivative(self.mesh,g.radial_points,g.degree,s) for s in (0,1,2)},
                               axial_derivative(self.mesh,g.axial_points,g.degree))
        dr,dz=self._derivatives
        nr=g.tensor(dr[0],density.cell_means,density.axial)
        nz=g.tensor(density.radial,density.cell_means,dz)
        fields=self.fields(density,p);v=fields['velocity'];pr=[];pz=[]
        for column,power in ((0,1),(2,2),(1,0)):
            rad=density.radial if power==0 else self.radial[power]
            pr.append(g.tensor(dr[power],p[...,column],density.axial))
            pz.append(g.tensor(rad,p[...,column],dz))
        pr=np.stack(pr,axis=-1);pz=np.stack(pz,axis=-1)
        pr[...,1]=pr[...,1]/g.r-fields['momentum'][...,1]/g.r
        pz[...,1]/=g.r
        vr=(pr-v*nr[...,None])/n[...,None];vz=(pz-v*nz[...,None])/n[...,None]
        curl=np.stack((-vz[...,1],vz[...,0]-vr[...,2],vr[...,1]+v[...,1]/g.r),axis=-1)
        div=vr[...,0]+v[...,0]/g.r+vz[...,2]
        return dict(radial=vr,axial=vz,curl=curl,divergence=div)
