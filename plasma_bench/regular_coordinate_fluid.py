"""Experimental joint Euler transport in x=r²/2 with dependent axis velocities."""
import numpy as np
from .lobatto_sbp import LobattoSBP
from .two_entropy_flux import TwoEntropyFlux
from .two_entropy_axial import TwoEntropyAxial
from .thermal_means import logarithmic_mean


class RegularCoordinateFluid:
    """U=(n,n*u,n*Omega,n*vz,n*se,n*si), u=vr/r; axis u/Omega slots inactive."""
    def __init__(self,radius,length,elements,axial_points,total_mass,*,degree=5,axial_order=6,gammas=(5/3,1.4)):
        self.radius=radius;self.length=length;self.elements=elements;self.nz=axial_points;self.degree=degree;self.mass=total_mass
        self.sbp=LobattoSBP(degree);self.thermal=TwoEntropyFlux(total_mass,gammas)
        if min(radius,length)<=0 or not np.isfinite([radius,length]).all() or elements<1 or axial_points<=axial_order:raise ValueError('Positive dimensions and resolved periodic axis required')
        self.coefficients=TwoEntropyAxial.coefficients[axial_order];self.dz=length/axial_points
        edges=.5*np.linspace(0,radius,elements+1)**2;self.dx=np.diff(edges)
        self.x=(edges[:-1,None]+.5*self.dx[:,None]*(1+self.sbp.nodes))[...,None];self.r=np.sqrt(2*self.x)
        self.z=(np.arange(axial_points)*self.dz-.5*length)[None,None,:]
        self.radial_weight=.5*self.dx[:,None,None]*self.sbp.weights[None,:,None]
        self.shape=(elements,degree+1,axial_points,6);self.volume=np.broadcast_to(2*np.pi*self.dz*self.radial_weight,self.shape[:-1])

    def fields(self,state):
        state=np.asarray(state,float)
        if state.shape!=self.shape or np.any(state[0,0,:,1:3]!=0):raise ValueError('Complete regular point state with inactive zero axis slots required')
        n=state[...,0]
        if np.any(n<=0) or not np.isfinite(state).all():raise ValueError('Positive finite regular state required')
        regular=state[...,1:4]/n[...,None];regular=regular.copy()
        regular[0,0,:,:2]=np.einsum('q,qzk->zk',self.sbp.axis_trace,regular[0,1:,:,:2])
        velocity=regular.copy();velocity[...,:2]*=self.r[...,None]
        physical=state.copy();physical[...,1:4]=n[...,None]*velocity;f=self.thermal.fields(physical)
        gradient=f['energy_gradient'].copy();gradient[...,1:3]*=self.r[...,None]
        return dict(f,regular_velocity=regular,energy_gradient=gradient,angular_momentum=2*self.x*n*regular[...,1])

    def pair(self,left,right,xleft,xright,radial):
        nl=logarithmic_mean(left['n'],right['n']);na=.5*(left['n']+right['n'])
        vl=left['regular_velocity'];vr=right['regular_velocity'];va=.5*(vl+vr)
        speed=xleft*vl[...,0]+xright*vr[...,0] if radial else va[...,2]
        gamma=nl*speed;pressure=.5*(np.sum(left['pressure'],axis=-1)+np.sum(right['pressure'],axis=-1))
        mom=gamma[...,None]*va;mom[...,0 if radial else 2]+=pressure/self.mass
        ha=.5*(left['specific_internal']+right['specific_internal']);hl=logarithmic_mean(left['specific_internal'],right['specific_internal'])
        sa=.5*(left['entropy']+right['entropy'])
        entropy=speed[...,None]*(nl[...,None]*sa+(self.thermal.gammas-1)*(na-nl)[...,None]+nl[...,None]*(ha/hl-1))
        return np.concatenate((gamma[...,None],mom,entropy),axis=-1)

    @staticmethod
    def subset(fields,index):
        return {k:fields[k][index] for k in ('n','regular_velocity','pressure','specific_internal','entropy')}

    def derivative(self,state):
        f=self.fields(state);q=self.sbp.q;left=self.subset(f,(slice(None),slice(None),None));right=self.subset(f,(slice(None),None,slice(None)))
        xi=self.x[:,:,None,:];xj=self.x[:,None,:,:];flux=self.pair(left,right,xi,xj,True)
        pressure=np.sum(f['pressure'],axis=-1);flux[...,1]-=.5*(pressure[:,:,None,:]+pressure[:,None,:,:])/self.mass
        residual=2*np.einsum('ij,eijzk->eizk',q,flux)
        residual[...,1]+=np.einsum('ij,eijz->eiz',q,pressure[:,None,:,:]-pressure[:,:,None,:])/self.mass
        u=f['regular_velocity'][...,0];omega=f['regular_velocity'][...,1]
        factor=(xj-xi)*logarithmic_mean(left['n'],right['n'])
        residual[...,1]-=np.einsum('ij,eijz->eiz',q,factor*(omega[:,:,None,:]**2-u[:,:,None,:]*u[:,None,:,:]))
        residual[...,2]+=2*u*np.einsum('ij,eijz->eiz',q,factor*.5*(omega[:,:,None,:]+omega[:,None,:,:]))
        boundary=self.subset(f,(slice(None),[0,-1]));xb=self.x[:,[0,-1]];self_flux=self.pair(boundary,boundary,xb,xb,True);star=np.zeros_like(self_flux)
        interface=self.pair(self.subset(f,(slice(None,-1),-1)),self.subset(f,(slice(1,None),0)),self.x[:-1,-1],self.x[1:,0],True)
        star[:-1,1]=interface;star[1:,0]=interface;star[0,0,:,1]=pressure[0,0]/self.mass;star[-1,1,:,1]=pressure[-1,-1]/self.mass
        residual[:,0]-=star[:,0]-self_flux[:,0];residual[:,-1]+=star[:,1]-self_flux[:,1]
        rate=-residual/self.radial_weight[...,None]
        center=self.subset(f,(slice(None),slice(None),slice(None)))
        for offset,coefficient in enumerate(self.coefficients,1):
            plus={k:np.roll(v,-offset,axis=2) for k,v in center.items()};minus={k:np.roll(v,offset,axis=2) for k,v in center.items()}
            rate-=2*coefficient/self.dz*(self.pair(center,plus,self.x,self.x,False)-self.pair(center,minus,self.x,self.x,False))
        rate[0,0,:,1:3]=0
        if not np.isfinite(rate).all():raise ValueError('Nonfinite regular coordinate rate')
        return rate,dict(physical_energy_rate=np.sum(f['energy_gradient']*rate,axis=-1))

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
        k1,_=self.derivative(state);k2,_=self.derivative(state+.5*dt*k1);k3,_=self.derivative(state+.5*dt*k2);k4,_=self.derivative(state+dt*k3)
        result=state+dt*(k1+2*k2+2*k3+k4)/6;self.fields(result);return result
