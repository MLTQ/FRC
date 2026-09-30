"""Cylindrical acoustic spectrum, original RHS controls, and RK4 time errors."""
import numpy as np
from scipy.linalg import eig,expm
from scipy.special import jn_zeros
from scipy.constants import elementary_charge as E
from plasma_bench.cylindrical_entropy_fluid import CylindricalEntropyFluid
from .log_fluid_reference import LogFluidReference
from .cylindrical_entropy_reference import rates


def calculate(spec):
    ref=LogFluidReference(gamma_i=5/3)
    d=CylindricalEntropyFluid(ref.radius,2*ref.half_length,spec['elements'],8,ref.mass,degree=spec['degree'],gammas=(5/3,5/3))
    count=d.elements*(d.degree+1);q=d.sbp.extended;r=np.concatenate((d.r,d.boundary_radius),axis=1)[...,0]
    area=.5*(r[:,:,None]+r[:,None,:]);rb=d.boundary_radius[...,0];sign=np.array([-1.,1.])
    def linear(x):
        fields=x.reshape(2,d.elements,d.degree+1);ext=np.concatenate((fields,np.einsum('bq,keq->keb',d.sbp.trace,fields)),axis=2)
        density,velocity=ext
        mass=np.einsum('ij,eij->ei',q,area*(velocity[:,:,None]+velocity[:,None,:]))
        force=np.einsum('ij,eij->ei',q,area*(density[:,None,:]-density[:,:,None]))
        boundary=ext[:,:,-2:];star=boundary.copy();star[1,0,0]=0;star[1,-1,1]=0
        interface=.5*(boundary[:,:-1,1]+boundary[:,1:,0]);star[:,:-1,1]=interface;star[:,1:,0]=interface
        residual=np.stack((mass,force));surface=sign*rb*(star[::-1]-boundary[::-1])
        lifted=np.einsum('qi,kei->keq',d.sbp.lift,residual)+np.einsum('bq,keb->keq',d.sbp.trace,surface)
        return (-ref.radius*lifted/d.radial_weight[...,0]).ravel()
    a=np.stack([linear(x) for x in np.eye(2*count)],axis=1)
    h=np.diag(np.tile(d.radial_weight[...,0].ravel(),2));values,vectors=eig(a)
    positive=np.flatnonzero(values.imag>1e-6);positive=positive[np.argsort(values[positive].imag)]
    analytic=jn_zeros(1,3);measured=values[positive[:3]].imag
    mode=vectors[:,positive[0]];x=mode.real if np.linalg.norm(mode.real)>.1 else mode.imag;x=x/np.max(abs(x))
    s=np.log(ref.pressure*np.array([.6,.4])/(ref.number_density*E*(ref.number_density/1e20)**(2/3)))
    def nonlinear(delta):
        u=np.zeros(d.shape);u[...,0]=ref.number_density*(1+delta*x[:count].reshape(d.elements,d.degree+1,1))
        u[...,1]=u[...,0]*ref.sound*delta*x[count:].reshape(d.elements,d.degree+1,1);u[...,4:]=u[...,0,None]*s
        du,_=d.derivative(u);f=rates(u,du,d.r,ref)
        return np.r_[f['n'].mean(axis=2).ravel()/ref.number_density,f['velocity'][...,0].mean(axis=2).ravel()/ref.sound]*ref.radius/ref.sound
    controls=[]
    for epsilon in (1e-3,1e-4,1e-5):
        actual=(nonlinear(epsilon)-nonlinear(-epsilon))/(2*epsilon)
        controls.append(dict(amplitude=epsilon,relative_error=float(np.linalg.norm(actual-a@x)/np.linalg.norm(a@x))))
    duration=4*np.pi/measured[0];exact=expm(a*duration)@x;e0=x@h@x;temporal=[]
    base=4*int(np.ceil(duration*max(abs(values))/.6/4))
    for steps in (base,2*base,4*base):
        y=x.copy();dt=duration/steps
        for _ in range(steps):
            k1=a@y;k2=a@(y+.5*dt*k1);k3=a@(y+.5*dt*k2);k4=a@(y+dt*k3);y+=dt*(k1+2*k2+2*k3+k4)/6
        temporal.append(dict(steps=steps,relative_error=float(np.linalg.norm(y-exact)/np.linalg.norm(exact)),energy_ratio=float(y@h@y/e0)))
    record=dict(spec=spec,maximum_real_scaled=float(max(values.real)),maximum_frequency_scaled=float(max(abs(values.imag))),
                energy_adjoint_defect=float(np.linalg.norm(h@a+a.T@h)/np.linalg.norm(h@a)),
                first_three_frequency_relative_errors=(measured/analytic-1).tolist(),nonlinear_difference_controls=controls,
                matrix_exponential_energy_ratio=float(exact@h@exact/e0),rk4_controls=temporal)
    return record,dict(operator=a,energy_metric=h,eigenvalues=values,initial=x,exact_final=exact,analytic_frequencies=analytic,measured_frequencies=measured)
