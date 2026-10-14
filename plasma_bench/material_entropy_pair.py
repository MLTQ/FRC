"""Species-relative entropy transport and its pressure work on material carriers."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME
from scipy.sparse import diags, kron
from scipy.sparse.linalg import splu
from .compatible_polynomials import bases


class MaterialEntropyPair:
    """Physical point entropies, fixed particle weights and compatible gradients.

    Species order is electron, ion. No collision or heat source is implied.
    The current geometry and drift must be from the same material state.
    """
    def __init__(self, fluid, space):
        self.fluid=fluid;self.space=space;self.particles=fluid.particles.ravel()
        mesh=space.mesh;q=space.order
        rr=space.reference_r.reshape(mesh.nr,mesh.nz,q,q)[:,0,:,0].ravel()
        zz=space.reference_z.reshape(mesh.nr,mesh.nz,q,q)[0,:,0,:].ravel()
        nr=bases(.5*mesh.rf**2,.5*rr**2,space.degree)[0]
        nz=bases(mesh.zf,zz,space.degree)[0]
        order=np.arange(len(rr)*len(zz)).reshape(mesh.nr,q,mesh.nz,q).transpose(0,2,1,3).ravel()
        self.scalar=kron(nr,nz,format='csr')[order]
        self.weighted=self.scalar.T@diags(self.particles)
        self.factor=splu((self.weighted@self.scalar).tocsc())
        self.fractions=(-fluid.ion_mass/fluid.mass,ME/fluid.mass)

    def project(self,value):
        value=np.asarray(value,float)
        if value.shape!=(len(self.particles),) or not np.isfinite(value).all():
            raise ValueError('Finite scalar on the actual material quadrature required')
        return self.factor.solve(self.weighted@value)

    def gradient(self,value,geometry):
        coefficient=self.space.gradient@self.project(value)
        return np.column_stack((geometry['er']@coefficient,geometry['ez']@coefficient))

    def gradient_transpose(self,vector,geometry):
        load=geometry['er'].T@vector[:,0]+geometry['ez'].T@vector[:,1]
        return self.weighted.T@self.factor.solve(self.space.gradient.T@load,trans='T')

    def transport(self,entropy,drift,geometry):
        """K*s with K=(A-A.T)/2 and A*s=dN*w.grad(P*s)."""
        forward=self.particles*np.sum(drift*self.gradient(entropy,geometry),axis=-1)
        reverse=self.gradient_transpose((self.particles*entropy)[:,None]*drift,geometry)
        return .5*(forward-reverse)

    def pair(self,state,geometry,drift):
        entropy=state.entropy.reshape(-1,2)
        internal=geometry['fields']['specific_internal'].reshape(-1,2)
        drift=np.asarray(drift,float)
        if drift.shape!=(len(self.particles),2) or not np.isfinite(drift).all():
            raise ValueError('Actual finite meridional carrier drift required')
        rate=np.empty_like(entropy);forces=[];potentials=[]
        for k,(fraction,gamma) in enumerate(zip(self.fractions,self.fluid.gammas)):
            s=entropy[:,k];h=internal[:,k];sa=float(s[0]);ha=float(h[0])
            ds=s-sa;dh=h-ha
            rate[:,k]=-fraction*self.transport(ds,drift,geometry)/self.particles
            reduced=-.5*fraction/E*(dh[:,None]*self.gradient(ds,geometry)-ds[:,None]*self.gradient(dh,geometry))
            forces.append(geometry['er'].T@(self.particles*reduced[:,0])+geometry['ez'].T@(self.particles*reduced[:,1]))
            phi=fraction/E*self.project(gamma*h-.5*s*h-.5*ha*s+.5*sa*h)
            potentials.append(self.space.gradient@phi)
        if not np.isfinite(rate).all():raise ValueError('Nonfinite material species transport; no clipping')
        return dict(entropy_rate=rate.reshape(state.entropy.shape),reduced_forcing=sum(forces),
            potential_electric=sum(potentials),
            thermal_power=float(np.sum(self.particles[:,None]*internal*rate)),
            entropy_rates_over_kB=np.sum(self.particles[:,None]*rate/(self.fluid.gammas-1),axis=0))

    def heating_rate(self,geometry,heating):
        fields=geometry['fields'];n=fields['n'].ravel();internal=fields['specific_internal'].reshape(-1,2)
        return ((heating/n[:,None])/internal).reshape(self.fluid.particles.shape+(2,))

    def potential_gradient(self,value):return self.space.gradient@self.project(value)
