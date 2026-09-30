"""Periodic axial flux-difference verification of joint two-entropy transport."""
import numpy as np
from .two_entropy_flux import TwoEntropyFlux


class TwoEntropyAxial:
    coefficients={2:(.5,),4:(2/3,-1/12),6:(3/4,-3/20,1/60)}
    def __init__(self,points,length,total_mass,*,gammas=(5/3,1.4),order=6):
        if order not in self.coefficients or points<=order or not np.isfinite(length) or length<=0:
            raise ValueError('Positive periodic length, sufficient points and order2/4/6 required')
        self.points=int(points);self.length=float(length);self.dx=length/points;self.order=order
        self.flux=TwoEntropyFlux(total_mass,gammas);self.z=np.arange(points)*self.dx-.5*length

    def validate(self,state):
        if np.shape(state)!=(self.points,6):raise ValueError('One six-component point state per periodic node required')
        return self.flux.fields(state)

    def derivative(self,state):
        physical=self.validate(state);rate=np.zeros_like(state);energy=np.zeros(self.points)
        for offset,coefficient in enumerate(self.coefficients[self.order],1):
            plus,gplus=self.flux.pair(state,np.roll(state,-offset,axis=0))
            minus,gminus=self.flux.pair(state,np.roll(state,offset,axis=0))
            rate-=2*coefficient/self.dx*(plus-minus);energy-=2*coefficient/self.dx*(gplus-gminus)
        actual=np.sum(physical['energy_gradient']*rate,axis=-1)
        return rate,dict(physical_energy_rate=actual,energy_flux_rate=energy)

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Finite positive step required')
        k1,_=self.derivative(state);k2,_=self.derivative(state+.5*dt*k1)
        k3,_=self.derivative(state+.5*dt*k2);k4,_=self.derivative(state+dt*k3)
        result=state+dt*(k1+2*k2+2*k3+k4)/6;self.validate(result);return result
