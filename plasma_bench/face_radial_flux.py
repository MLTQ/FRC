"""Conservative axial/angular fluxes with the actual radial face momentum."""
import numpy as np
from scipy.constants import elementary_charge as E
from .log_fluid_flux import LogFluidFlux


class FaceRadialFlux(LogFluidFlux):
    def __init__(self,geometry,total_mass,kinematics):
        super().__init__(geometry,total_mass);self.kinematics=kinematics

    def traces(self,density,radial,momentum,species,states,normal):
        m=self.mesh;g=self.geometry
        if normal=='radial':
            maps=self.rfaces;axial=density.axial;r=m.rf[:,None]
            q=self.kinematics.full_faces(radial)
            mr=np.divide(q,r,out=np.zeros_like(q),where=r>0)@axial.T
            entropies=[np.asarray(s@self.znodes.T) for s in states]
        elif normal=='axial':
            maps=self.radial;axial=self.zfaces;r=g.radial_points[:,None]
            mr=np.asarray(self.kinematics.radial@radial@axial.T)
            entropies=[np.asarray(self.rnodes@s) for s in states]
        else:raise ValueError('Radial or axial normal required')
        n=np.asarray(maps[0]@density.cell_means@axial.T)
        if not np.isfinite(n).all() or np.any(n<=0):raise ValueError('Nonpositive fluid face density; no floor')
        mz=np.asarray(maps[0]@momentum[...,0]@axial.T)
        angular=np.asarray(maps[2]@momentum[...,1]@axial.T)
        toroidal=np.divide(angular,r,out=np.zeros_like(angular),where=r>0)
        vr,vz,vt=mr/n,mz/n,toroidal/n;pressures=[]
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            for thermal,s in zip(species,entropies):
                p=n*E*np.exp(s+(thermal.gamma-1)*np.log(n/thermal.density_scale))
                if not np.isfinite(p).all() or np.any(p<=0):raise ValueError('Nonpositive/nonfinite fluid face pressure')
                pressures.append(p)
        speed=vr if normal=='radial' else vz
        specific=np.stack((mz/n,angular/n),axis=-1)
        momentum_flux=self.mass*np.stack((mz*speed,angular*speed),axis=-1)
        if normal=='axial':momentum_flux[...,0]+=sum(pressures)
        kinetic=.5*self.mass*n*(vr*vr+vz*vz+vt*vt)
        enthalpy=sum(t.gamma/(t.gamma-1)*p for t,p in zip(species,pressures))
        energy_flux=(kinetic+enthalpy)*speed
        if not np.isfinite(momentum_flux).all() or not np.isfinite(energy_flux).all():raise ValueError('Fluid face flux exceeds finite range')
        return dict(momentum=momentum_flux,energy=energy_flux,particles=n*speed,specific=specific)

    def rate(self,model,radial,momentum,species,states):
        r=self.traces(model.density_map,radial,momentum,species,states,'radial')
        z=self.traces(model.density_map,radial,momentum,species,states,'axial')
        rm,re,rn=[self.integrate(r[k],'radial') for k in ('momentum','energy','particles')]
        zm,ze,zn=[self.integrate(z[k],'axial') for k in ('momentum','energy','particles')]
        return dict(momentum=-self.divergence(rm,zm)/self.mass,energy_account=-self.divergence(re,ze),
                    radial_momentum_flux=rm,axial_momentum_flux=zm,radial_particles=rn,axial_particles=zn,
                    radial_energy_flux=re,axial_energy_flux=ze,radial_specific_momentum=r['specific'])
