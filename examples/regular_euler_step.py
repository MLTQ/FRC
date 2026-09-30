"""Verification source work for the regular-coordinate physical state."""
import numpy as np
from .regular_euler_reference import energy


def source(driver,reference,time):
    value=reference.regular_source(driver.x,driver.z,time);value[0,0,:,1:3]=0.;return value


def integrals(driver,value):
    weighted=value.copy();weighted[...,2]*=2*driver.x
    return np.sum(driver.volume[...,None]*weighted,axis=(0,1,2))


def step(driver,state,reference,time,dt,manufactured):
    if not manufactured:return driver.step(state,dt),0.,np.zeros(6)
    if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
    def rhs(u,t):
        rate,_=driver.derivative(u);forcing=source(driver,reference,t)
        work=float(np.sum(driver.volume*np.imag(energy(u+1j*1e-28*forcing,driver.x,reference)))/1e-28)
        return rate+forcing,work,integrals(driver,forcing)
    k1,w1,s1=rhs(state,time);k2,w2,s2=rhs(state+.5*dt*k1,time+.5*dt)
    k3,w3,s3=rhs(state+.5*dt*k2,time+.5*dt);k4,w4,s4=rhs(state+dt*k3,time+dt)
    result=state+dt*(k1+2*k2+2*k3+k4)/6;driver.fields(result)
    return result,dt*(w1+2*w2+2*w3+w4)/6,dt*(s1+2*s2+2*s3+s4)/6
