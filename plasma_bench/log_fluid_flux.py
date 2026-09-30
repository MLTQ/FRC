"""Shared finite-volume momentum and enthalpy fluxes for the physical log state."""
import numpy as np
from scipy.constants import elementary_charge as E
from .compatible_polynomials import bases
from .compatible_velocity import radial_moments
from .normal_momentum_faces import radial_faces,axial_faces


class LogFluidFlux:
    def __init__(self,geometry,total_mass):
        self.geometry=geometry;self.mesh=geometry.mesh;self.mass=total_mass
        g=geometry;m=self.mesh;p=g.degree
        self.radial={k:radial_moments(m,g.radial_points,p,power=k) for k in (0,1,2)}
        self.rfaces={k:radial_faces(m,p,power=k) for k in (0,1,2)}
        self.zfaces=axial_faces(m,p)
        self.rnodes=bases(.5*m.rf**2,.5*g.radial_points**2,p)[0]
        self.znodes=bases(m.zf,g.axial_points,p,periodic=True)[0]
        _,self.gauss_weights=np.polynomial.legendre.leggauss(g.order)

    def traces(self,density,momentum,species,states,normal):
        m=self.mesh;g=self.geometry;degree=g.degree
        if normal=='radial':
            maps=self.rfaces;axial=density.axial;r=m.rf[:,None]
            n=np.asarray(maps[0]@density.cell_means@axial.T)
            entropies=[np.asarray(s@self.znodes.T) for s in states]
        elif normal=='axial':
            maps=self.radial;axial=self.zfaces;r=g.radial_points[:,None]
            n=np.asarray(maps[0]@density.cell_means@axial.T)
            entropies=[np.asarray(self.rnodes@s) for s in states]
        else:raise ValueError('Radial or axial flux normal required')
        if not np.isfinite(n).all() or np.any(n<=0):raise ValueError('Nonpositive fluid face density; no floor')
        mr=np.asarray(maps[1]@momentum[...,0]@axial.T)
        mz=np.asarray(maps[0]@momentum[...,1]@axial.T)
        angular=np.asarray(maps[2]@momentum[...,2]@axial.T)
        if normal=='radial':mr[[0,-1]]=0.
        toroidal=np.divide(angular,r,out=np.zeros_like(angular),where=r>0)
        vr,vz,vt=mr/n,mz/n,toroidal/n
        pressures=[]
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            for thermal,s in zip(species,entropies):
                pressure=n*E*np.exp(s+(thermal.gamma-1)*np.log(n/thermal.density_scale))
                if not np.isfinite(pressure).all() or np.any(pressure<=0):raise ValueError('Nonpositive/nonfinite fluid face pressure')
                pressures.append(pressure)
        p=sum(pressures);speed=vr if normal=='radial' else vz
        momentum_flux=self.mass*np.stack((mr*speed,mz*speed,angular*speed),axis=-1)
        momentum_flux[...,0 if normal=='radial' else 1]+=p
        kinetic=.5*self.mass*n*(vr*vr+vz*vz+vt*vt)
        enthalpy=sum(t.gamma/(t.gamma-1)*ps for t,ps in zip(species,pressures))
        energy_flux=(kinetic+enthalpy)*speed
        if not np.isfinite(momentum_flux).all() or not np.isfinite(energy_flux).all():raise ValueError('Fluid face flux exceeds finite range')
        return momentum_flux,energy_flux,n*speed

    def integrate(self,values,normal):
        m=self.mesh;g=self.geometry;tail=values.shape[2:]
        if normal=='radial':
            v=values.reshape((m.nr+1,m.nz,g.order)+tail)
            w=(2*np.pi*m.rf[:,None,None]*.5*m.dz*self.gauss_weights[None,None,:])
            return np.sum(v*w[(...,)+(None,)*len(tail)],axis=2)
        v=values.reshape((m.nr,g.order,m.nz)+tail)
        w=g.wr.reshape(m.nr,g.order,1)
        return np.sum(v*w[(...,)+(None,)*len(tail)],axis=1)

    def divergence(self,radial,axial):
        diff=np.diff(radial,axis=0)+(np.roll(axial,-1,axis=1)-axial)
        return diff/self.mesh.volumes[(...,)+(None,)*(diff.ndim-2)]

    def rate(self,model,kinematics,momentum,species,states):
        rm,re,rn=self.traces(model.density_map,momentum,species,states,'radial')
        zm,ze,zn=self.traces(model.density_map,momentum,species,states,'axial')
        rm,re,rn=[self.integrate(v,'radial') for v in (rm,re,rn)]
        zm,ze,zn=[self.integrate(v,'axial') for v in (zm,ze,zn)]
        force=-self.divergence(rm,zm)
        v=kinematics.fields(model.density_map,momentum)['velocity'];g=self.geometry
        p=sum(t.fields(s)['pressure'] for t,s in zip(species,states))
        hoop=model.average((p+self.mass*model.density*v[...,1]**2)/g.r)
        force[...,0]+=hoop
        return dict(momentum=force/self.mass,energy_account=-self.divergence(re,ze),
            radial_momentum_flux=rm,axial_momentum_flux=zm,radial_energy_flux=re,axial_energy_flux=ze,
            radial_particles=rn,axial_particles=zn,hoop_force=hoop)
