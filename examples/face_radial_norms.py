"""Independent common-flux radial momentum and physical energy sampling."""
import numpy as np
from scipy.constants import elementary_charge as E
from .log_fluid_norms import FluidSampler
from .bulk_thermal_norms import flux_basis


class FaceRadialSampler(FluidSampler):
    def __init__(self,mesh,degree=5,order=14):
        super().__init__(mesh,degree,order)
        self.Q=flux_basis(.5*mesh.rf**2,.5*self.r[:,0]**2,degree)[:,1:-1]/self.r

    def momentum(self,radial,values):
        return np.stack((self.Q@radial@self.Rz.T,
                         (self.regular[2]@values[...,1]@self.Rz.T)/self.r,
                         self.density(values[...,0])),axis=-1)

    def physical(self,state,reference):
        n=self.density(state.density);p=self.momentum(state.radial,state.momentum);v=p/n[...,None]
        f=dict(n=n,velocity=v,kinetic=.5*reference.mass*np.sum(p*v,axis=-1))
        for label,short,gamma in (('electron','e',reference.gamma_e),('ion','i',reference.gamma_i)):
            s=self.scalar(getattr(state,label));temp=np.exp(s)*(n/1e20)**(gamma-1)
            f[label]=s;f['T'+short]=temp;f['p'+short]=n*E*temp
        f['energy']=f['kinetic']+f['pe']/(reference.gamma_e-1)+f['pi']/(reference.gamma_i-1)
        return f

    def rates(self,state,rate,reference,time):
        f=self.physical(state,reference);_,_,dt=reference.derivatives(self.r,self.z,time)
        nt=self.density(rate.density);pt=self.momentum(rate.radial,rate.momentum)
        vt=(pt-f['velocity']*nt[...,None])/f['n'][...,None]
        energy=reference.mass*(np.sum(f['velocity']*pt,axis=-1)-.5*np.sum(f['velocity']**2,axis=-1)*nt)
        result=dict(n=self.metric(nt,dt['n']),velocity=self.metric(vt,dt['velocity']))
        result['velocity_components']={key:self.metric(vt[...,k],dt['velocity'][...,k])
                                       for k,key in enumerate(('radial','azimuthal','axial'))}
        for label,short,gamma in (('electron','e',reference.gamma_e),('ion','i',reference.gamma_i)):
            st=self.scalar(getattr(rate,label));pressure=f['p'+short]*(st+gamma*nt/f['n'])
            temperature=f['T'+short]*(st+(gamma-1)*nt/f['n'])
            result['p'+short]=self.metric(pressure,dt['p'+short]);result['T'+short]=self.metric(temperature,dt['T'+short])
            energy+=pressure/(gamma-1)
        result['energy']=self.metric(energy,dt['energy'])
        return result,self.cells(energy)
