"""Verification-only manufactured forcing with actual physical source work."""
import numpy as np
from .cylindrical_entropy_reference import rates


def step(driver,state,reference,time,dt,manufactured):
    if not manufactured:return driver.step(state,dt),0.,np.zeros(6)
    if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite verification step required')
    def rhs(u,t):
        rate,_=driver.derivative(u);source=reference.source(driver.r,driver.z,t)
        # The energy source is evaluated at the actual numerical state,
        # using the original kinetic/EOS chain rule, not the exact target H.
        work=float(np.sum(driver.volume*rates(u,source,driver.r,reference)['energy']))
        totals=np.sum(driver.volume[...,None]*source,axis=(0,1,2))
        return rate+source,work,totals
    k1,w1,s1=rhs(state,time);k2,w2,s2=rhs(state+.5*dt*k1,time+.5*dt)
    k3,w3,s3=rhs(state+.5*dt*k2,time+.5*dt);k4,w4,s4=rhs(state+dt*k3,time+dt)
    result=state+dt*(k1+2*k2+2*k3+k4)/6;driver.extended(result)
    return result,dt*(w1+2*w2+2*w3+w4)/6,dt*(s1+2*s2+2*s3+s4)/6
