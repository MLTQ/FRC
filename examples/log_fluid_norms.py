"""Independent physical reconstruction and norms for the evolving fluid state."""
import numpy as np
from scipy.constants import elementary_charge as E
from .physical_scalar_sampler import PhysicalScalarSampler


def regular_map(mesh,points,degree,power):
    """Independent Legendre-moment fit, using higher-order physical integration."""
    out=np.zeros((len(points),mesh.nr));q,w=np.polynomial.legendre.leggauss(20)
    for row,r in enumerate(points):
        cell=min(np.searchsorted(mesh.rf,r,side='right')-1,mesh.nr-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            first=min(max(cell-offset,0),mesh.nr-degree-1);ids=np.arange(first,first+degree+1)
            lo,hi=mesh.rf[first]**2,mesh.rf[first+degree+1]**2
            rq=mesh.r[ids,None]+.5*mesh.dr*q;eta=(2*rq*rq-lo-hi)/(hi-lo)
            pol=np.polynomial.legendre.legvander(eta,degree)
            weights=.5*mesh.dr*rq*w/mesh.radial_measure[ids,None]
            moments=np.sum(weights[...,None]*(rq/mesh.radius)[...,None]**power*pol,axis=1)
            target=(r/mesh.radius)**power*np.polynomial.legendre.legvander((2*r*r-lo-hi)/(hi-lo),degree)[0]
            out[row,ids]+=.5*np.linalg.solve(moments.T,target)
    return out


class FluidSampler(PhysicalScalarSampler):
    def __init__(self,mesh,degree=5,order=14):
        super().__init__(mesh,degree,order)
        self.regular={p:regular_map(mesh,self.r[:,0],degree,p) for p in (1,2)}

    def momentum(self,values):
        return np.stack((self.regular[1]@values[...,0]@self.Rz.T,
                         (self.regular[2]@values[...,2]@self.Rz.T)/self.r,
                         self.density(values[...,1])),axis=-1)

    def physical(self,state,reference):
        n=self.density(state.density);p=self.momentum(state.momentum);v=p/n[...,None]
        f=dict(n=n,velocity=v,kinetic=.5*reference.mass*np.sum(p*v,axis=-1))
        for label,short,gamma in (('electron','e',reference.gamma_e),('ion','i',reference.gamma_i)):
            s=self.scalar(getattr(state,label));temp=np.exp(s)*(n/1e20)**(gamma-1)
            f[label]=s;f['T'+short]=temp;f['p'+short]=n*E*temp
        f['energy']=f['kinetic']+f['pe']/(reference.gamma_e-1)+f['pi']/(reference.gamma_i-1)
        return f

    def measure(self,state,reference,time):
        f=self.physical(state,reference);exact=reference.fields(self.r,self.z,time)
        result={k:self.metric(f[k],exact[k]) for k in ('n','velocity','pe','pi','Te','Ti','energy')}
        result['density_signal']=self.metric(f['n']-reference.number_density,exact['n']-reference.number_density)
        result['energy_integral']=float(np.sum(self.weight*f['energy']))
        result['entropy_integrals']=[float(np.sum(self.weight*f['n']*f[k])) for k in ('electron','ion')]
        return result,self.cells(f['energy'])

    def rates(self,state,rate,reference,time):
        f=self.physical(state,reference);_,_,dt=reference.derivatives(self.r,self.z,time)
        nt=self.density(rate.density);pt=self.momentum(rate.momentum)
        vt=(pt-f['velocity']*nt[...,None])/f['n'][...,None]
        kinetic=reference.mass*(np.sum(f['velocity']*pt,axis=-1)-.5*np.sum(f['velocity']**2,axis=-1)*nt)
        result=dict(n=self.metric(nt,dt['n']),velocity=self.metric(vt,dt['velocity']))
        energy=kinetic.copy()
        for label,short,gamma in (('electron','e',reference.gamma_e),('ion','i',reference.gamma_i)):
            st=self.scalar(getattr(rate,label));pressure=f['p'+short]*(st+gamma*nt/f['n'])
            temperature=f['T'+short]*(st+(gamma-1)*nt/f['n'])
            result['p'+short]=self.metric(pressure,dt['p'+short]);result['T'+short]=self.metric(temperature,dt['T'+short])
            energy+=pressure/(gamma-1)
        result['energy']=self.metric(energy,dt['energy'])
        return result,self.cells(energy)
