"""Cylindrical joint transport with conservative angular momentum and geometric work."""
import numpy as np
from .gauss_sbp import GaussSBP
from .cylindrical_entropy_state import CylindricalEntropyState
from .two_entropy_axial import TwoEntropyAxial
from .thermal_means import logarithmic_mean


class CylindricalEntropyFluid:
    def __init__(self,radius,length,elements,axial_points,total_mass,*,degree=5,axial_order=6,gammas=(5/3,1.4)):
        if min(radius,length)<=0 or not np.isfinite([radius,length,total_mass]).all() or elements<1:
            raise ValueError('Positive finite cylindrical dimensions and element count required')
        if axial_order not in TwoEntropyAxial.coefficients or axial_points<=axial_order:
            raise ValueError('Enough periodic axial points for order2/4/6 required')
        self.radius=radius;self.length=length;self.elements=elements;self.nz=axial_points;self.degree=degree
        self.sbp=GaussSBP(degree);self.mass=total_mass;self.state_map=CylindricalEntropyState(total_mass,gammas)
        self.dr=radius/elements;self.dz=length/axial_points;self.coefficients=TwoEntropyAxial.coefficients[axial_order]
        edges=np.linspace(0,radius,elements+1);self.boundary_radius=np.stack((edges[:-1],edges[1:]),axis=1)[...,None]
        self.r=(edges[:-1,None]+.5*self.dr*(1+self.sbp.nodes))[...,None]
        self.z=(np.arange(axial_points)*self.dz-.5*length)[None,None,:]
        self.radial_weight=.5*self.dr*self.sbp.weights[None,:,None]*self.r
        self.shape=(elements,degree+1,axial_points,6)
        self.volume=np.broadcast_to(2*np.pi*self.dz*self.radial_weight,self.shape[:-1])

    def fields(self,state):
        if np.shape(state)!=self.shape:raise ValueError('Cylindrical point state has the wrong shape')
        return self.state_map.fields(state,self.r)

    def extended(self,state):
        volume=self.fields(state);w=np.einsum('bq,eqzk->ebzk',self.sbp.trace,volume['energy_gradient'])
        boundary,bf=self.state_map.from_energy_variables(w,self.boundary_radius)
        physical=np.concatenate((volume['physical_state'],boundary),axis=1)
        radius=np.concatenate((self.r,self.boundary_radius),axis=1)
        omega=np.concatenate((volume['omega'],bf['omega']),axis=1)
        return volume,boundary,physical,radius,omega

    def radial_rate(self,state):
        volume,boundary,u,r,omega=self.extended(state);q=self.sbp.extended
        left=u[:,:,None,:,:];right=u[:,None,:,:,:];rl=r[:,:,None,:];rr=r[:,None,:,:]
        area=.5*(rl+rr);flux=self.state_map.pair(left,right,rl,rr,0)
        f=self.state_map.thermal.fields(u);pressure=np.sum(f['pressure'],axis=-1)
        # Remove pressure before differencing; constant physical pressure has
        # exactly zero gradient, without subtracting large geometric sources.
        flux[...,1]-=.5*(pressure[:,:,None,:]+pressure[:,None,:,:])/self.mass
        residual=2*np.einsum('ij,eijzk->eizk',q,area[...,None]*flux)
        jump=pressure[:,None,:,:]-pressure[:,:,None,:]
        residual[...,1]+=np.einsum('ij,eijz->eiz',q,area*jump)/self.mass
        nlog=logarithmic_mean(f['n'][:,:,None,:],f['n'][:,None,:,:])
        transverse=.5*(f['velocity'][:,:,None,:,1]+f['velocity'][:,None,:,:,1])
        angular=.5*(omega[:,:,None,:]+omega[:,None,:,:])
        geometric=area*nlog*(rr-rl)*transverse*angular
        residual[...,1]-=np.einsum('ij,eijz->eiz',q,geometric)
        rb=self.boundary_radius;self_flux=self.state_map.pair(boundary,boundary,rb,rb,0)
        star=np.zeros_like(self_flux)
        interface=self.state_map.pair(boundary[:-1,1],boundary[1:,0],rb[:-1,1],rb[1:,0],0)
        star[:-1,1]=interface;star[1:,0]=interface
        for element,side in ((0,0),(-1,1)):
            wall=self.state_map.thermal.fields(boundary[element,side])
            star[element,side,:,1]=np.sum(wall['pressure'],axis=-1)/self.mass
        surface=np.array([-1.,1.])[None,:,None,None]*rb[...,None]*(star-self_flux)
        integrated=np.einsum('qi,eizk->eqzk',self.sbp.lift,residual)+np.einsum('bq,ebzk->eqzk',self.sbp.trace,surface)
        return -integrated/self.radial_weight[...,None],volume

    def derivative(self,state):
        rate,volume=self.radial_rate(state);physical=volume['physical_state']
        for offset,coefficient in enumerate(self.coefficients,1):
            plus=self.state_map.pair(physical,np.roll(physical,-offset,axis=2),self.r,self.r,2)
            minus=self.state_map.pair(physical,np.roll(physical,offset,axis=2),self.r,self.r,2)
            rate-=2*coefficient/self.dz*(plus-minus)
        if not np.isfinite(rate).all():raise ValueError('Cylindrical physical rate exceeds finite range')
        energy=np.sum(volume['energy_gradient']*rate,axis=-1)
        return rate,dict(physical_energy_rate=energy)

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Finite positive step required')
        k1,_=self.derivative(state);k2,_=self.derivative(state+.5*dt*k1)
        k3,_=self.derivative(state+.5*dt*k2);k4,_=self.derivative(state+dt*k3)
        result=state+dt*(k1+2*k2+2*k3+k4)/6;self.extended(result);return result
