"""Regular two-species fluid in a finite cylinder with impermeable end caps."""
import numpy as np
from .lobatto_sbp import LobattoSBP
from .strong_wall_plasma import StrongWallFluid


class FiniteCylinderFluid(StrongWallFluid):
    def __init__(self, radius, length, elements, axial_points, total_mass, *, degree=5, gammas=(5/3,5/3)):
        super().__init__(radius,length,elements,axial_points,total_mass,degree=degree,axial_order=2,gammas=gammas)
        self.axial_sbp=LobattoSBP(axial_points-1)
        self.z=(.5*length*self.axial_sbp.nodes)[None,None,:]
        self.axial_weight=.5*length*self.axial_sbp.weights
        self.volume=np.broadcast_to(2*np.pi*self.radial_weight*self.axial_weight[None,None,:],self.shape[:-1])
        # Reuse the original radial transport; there is no periodic z term.
        self.coefficients=()

    def fields(self, state):
        result=super().fields(state)
        if np.any(np.asarray(state)[...,[0,-1],3] != 0):
            raise ValueError('Inactive zero end-normal momentum slots required')
        return result

    def derivative(self, state):
        rate,_=super().derivative(state); f=self.fields(state)
        left=self.subset(f,(slice(None),slice(None),slice(None),None))
        right=self.subset(f,(slice(None),slice(None),None,slice(None)))
        x=self.x[...,None]
        flux=self.pair(left,right,x,x,False)
        pressure=np.sum(f['pressure'],axis=-1)
        # Difference-first pressure derivative avoids differentiating its large
        # uniform background through an imperfect floating row sum.
        flux[...,3]-=.5*(pressure[..., :,None]+pressure[...,None,:])/self.mass
        residual=2*np.einsum('ij,erijk->erik',self.axial_sbp.q,flux)
        residual[...,3]+=np.einsum('ij,erij->eri',self.axial_sbp.q,
            pressure[...,None,:]-pressure[..., :,None])/self.mass
        boundary=self.subset(f,(slice(None),slice(None),[0,-1]))
        self_flux=self.pair(boundary,boundary,self.x,self.x,False)
        star=np.zeros_like(self_flux); star[...,3]=pressure[...,[0,-1]]/self.mass
        residual[...,0,:]-=star[...,0,:]-self_flux[...,0,:]
        residual[...,-1,:]+=star[...,-1,:]-self_flux[...,-1,:]
        rate-=residual/self.axial_weight[None,None,:,None]
        # Normal reactions enforce the declared physical wall coordinates and
        # do zero work because the corresponding velocity is exactly zero.
        reaction=-rate[...,[0,-1],3].copy(); rate[...,[0,-1],3]=0.
        rate[-1,-1,:,1]=0.; rate[0,0,:,1:3]=0.
        if not np.isfinite(rate).all():raise ValueError('Nonfinite finite-cylinder fluid rate')
        area=2*np.pi*self.radial_weight[...,0]
        wall_force=float(np.sum(area*(pressure[...,0]-pressure[...,-1])))
        reaction_force=float(self.mass*np.sum(self.volume[...,[0,-1]]*reaction))
        return rate,dict(physical_energy_rate=np.sum(f['energy_gradient']*rate,axis=-1),
            end_pressure_force_n=wall_force,end_constraint_force_n=reaction_force,
            axial_momentum_rate_n=float(self.mass*np.sum(self.volume*rate[...,3])))
