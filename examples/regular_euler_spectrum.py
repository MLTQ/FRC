"""Original radial acoustic linearization for the active regular-coordinate DOFs."""
import numpy as np
from scipy.special import jn_zeros
from scipy.linalg import eig,expm
from scipy.constants import elementary_charge as E
from .regular_euler_reference import RegularEulerReference
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid


def calculate(elements,degree):
    ref=RegularEulerReference(gamma_i=5/3);d=RegularCoordinateFluid(ref.radius,2*ref.half_length,elements,8,ref.mass,degree=degree,gammas=(5/3,5/3))
    x=d.x[...,0];size=x.size;q=d.sbp.q;jac=d.radial_weight[...,0];sign=np.array([-1.,1.])
    def action(y):
        density=y[:size].reshape(x.shape);speed=np.zeros(size);speed[1:]=y[size:];speed=speed.reshape(x.shape)
        speed[0,0]=d.sbp.axis_trace@speed[0,1:]
        bd=density[:,[0,-1]];bs=speed[:,[0,-1]];sd=bd.copy();ss=bs.copy();ss[0,0]=0.;ss[-1,1]=0.
        sd[:-1,1]=sd[1:,0]=.5*(bd[:-1,1]+bd[1:,0]);ss[:-1,1]=ss[1:,0]=.5*(bs[:-1,1]+bs[1:,0])
        dn=-2*np.einsum('ij,ej->ei',q,x*speed);da=-ref.radius**2*np.einsum('ij,ej->ei',q,density)
        dn[:,[0,-1]]-=sign*2*x[:,[0,-1]]*(ss-bs);da[:,[0,-1]]-=ref.radius**2*sign*(sd-bd)
        return np.r_[(dn/jac).ravel(),(da/jac).ravel()[1:]]
    a=np.stack([action(y) for y in np.eye(2*size-1)],axis=1)
    weights=np.r_[jac.ravel(),(2*x/ref.radius**2*jac).ravel()[1:]];h=np.diag(weights);values,vectors=eig(a)
    indices=np.flatnonzero(values.imag>1e-6);indices=indices[np.argsort(values[indices].imag)];positive=values[indices].imag;analytic=jn_zeros(1,3)
    mode=vectors[:,indices[0]];y=mode.real+mode.imag;y=y/np.max(abs(y))
    entropy=np.log(ref.pressure*np.array([.6,.4])/(ref.number_density*E*(ref.number_density/1e20)**(2/3)))
    def nonlinear(epsilon):
        state=np.zeros(d.shape);state[...,0]=ref.number_density*(1+epsilon*y[:size].reshape(x.shape))[...,None]
        speed=np.zeros(size);speed[1:]=epsilon*y[size:];speed=speed.reshape(x.shape)[...,None]
        state[...,1]=state[...,0]*ref.sound/ref.radius*speed;state[...,4:]=state[...,0,None]*entropy
        rate,_=d.derivative(state);ut=(rate[...,1]-state[...,1]*rate[...,0]/state[...,0])/state[...,0]
        return np.r_[np.mean(rate[...,0],axis=2).ravel()/ref.number_density*ref.radius/ref.sound,
                     np.mean(ut,axis=2).ravel()[1:]*ref.radius**2/ref.sound**2]
    controls=[]
    for epsilon in (1e-3,1e-4,1e-5):
        actual=(nonlinear(epsilon)-nonlinear(-epsilon))/(2*epsilon)
        controls.append(dict(amplitude=epsilon,relative_error=float(np.linalg.norm(actual-a@y)/np.linalg.norm(a@y))))
    duration=4*np.pi/positive[0];exact=expm(a*duration)@y;initial_energy=y@h@y;temporal=[]
    base=4*int(np.ceil(duration*np.max(abs(values))/.6/4))
    for count in (base,2*base,4*base):
        z=y.copy();dt=duration/count
        for _ in range(count):
            k1=a@z;k2=a@(z+.5*dt*k1);k3=a@(z+.5*dt*k2);k4=a@(z+dt*k3);z+=dt*(k1+2*k2+2*k3+k4)/6
        temporal.append(dict(steps=count,relative_error=float(np.linalg.norm(z-exact)/np.linalg.norm(exact)),energy_ratio=float(z@h@z/initial_energy)))
    return dict(elements=elements,degree=degree,maximum_real=float(np.max(values.real)),energy_adjoint_defect=float(np.linalg.norm(h@a+a.T@h)/np.linalg.norm(h@a)),
        frequency_errors=(positive[:3]/analytic-1).tolist(),nonlinear_controls=controls,rk4_controls=temporal),dict(operator=a,metric=h,eigenvalues=values,initial=y,exact_final=exact)
