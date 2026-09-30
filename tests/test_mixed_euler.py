import numpy as np
import pytest
from examples.mixed_euler_reference import MixedEulerReference
from examples.mixed_euler_step import step
from examples.mixed_euler_sampling import sample,grid,differences
from examples.cylindrical_entropy_reference import physical,rates
from plasma_bench.cylindrical_entropy_fluid import CylindricalEntropyFluid


def test_independent_material_and_conservative_manufactures_agree_with_wall():
    ref=MixedEulerReference();r=np.linspace(.001,1,47)[:,None]*ref.radius;z=np.linspace(-1,1,32)[None,:]*ref.half_length
    for t in (0.,.43*ref.duration,ref.duration):
        a=ref.source(r,z,t);b=ref.primitive_source(r,z,t)
        assert np.max(np.linalg.norm(a-b,axis=(0,1))/np.linalg.norm(a,axis=(0,1)))<2e-13
        wall=ref.fields(np.full_like(r,ref.radius),z,t);assert np.max(abs(wall['velocity'][...,0]))==0
        h=1e-28;f=ref.fields(np.full_like(r,ref.radius)+1j*h,z,t)
        gradient=np.imag(f['pe']+f['pi'])/h
        centrifugal=ref.mass*wall['n']*wall['velocity'][...,1]**2/ref.radius
        assert np.linalg.norm(gradient-centrifugal)/np.linalg.norm(centrifugal)<1e-14
        # Physical forcing work is the EOS derivative at the actual state.
        u=ref.state(r,z,t);wanted=np.imag(physical(u+1j*h*a,r,ref)['energy'])/h
        assert np.linalg.norm(wanted-rates(u,a,r,ref)['energy'])/np.linalg.norm(wanted)<2e-14


def test_common_grid_sampler_reproduces_regular_polynomial_and_fourier_fields():
    ref=MixedEulerReference();d=CylindricalEntropyFluid(ref.radius,2*ref.half_length,3,16,ref.mass,degree=3)
    r,z=np.broadcast_arrays(d.r,d.z);n=ref.number_density*(1+.2*(r/ref.radius)**2)*(1+.1*np.sin(np.pi*z/ref.half_length))
    u=np.zeros(d.shape);u[...,0]=n;u[...,1]=n*.1*ref.sound*r/ref.radius;u[...,2]=n*.2*ref.sound*r*r/ref.radius;u[...,4:]=n[...,None]
    rr,zz,volume=grid(ref.radius,2*ref.half_length,7,48,12);actual=sample(d,u,rr,len(zz));exact=ref.number_density*(1+.2*(rr[:,None]/ref.radius)**2)*(1+.1*np.sin(np.pi*zz[None,:]/ref.half_length))
    assert np.max(abs(actual[...,0]/exact-1))<5e-14
    assert np.max(abs(actual[...,2]/(exact*rr[:,None]**2)-.2*ref.sound/ref.radius))/(ref.sound/ref.radius)<5e-14


def test_time_dependent_manufacture_refines_and_source_step_rejects_atomically():
    ref=MixedEulerReference();rr,zz,volume=grid(ref.radius,2*ref.half_length,8,128,16);duration=.2*ref.duration;errors=[]
    for ne in (2,4):
        d=CylindricalEntropyFluid(ref.radius,2*ref.half_length,ne,8*ne,ref.mass);u=ref.state(d.r,d.z)
        count=80*ne;dt=duration/count
        for i in range(count):u,_,_=step(d,u,ref,i*dt,dt,True)
        metrics=differences(sample(d,u,rr,len(zz)),ref.state(rr[:,None],zz[None,:],duration),rr,volume,ref,ref.radius/8)
        errors.append(metrics['n']['relative_L2'])
    assert errors[0]/errors[1]>20
    before=u.copy()
    with pytest.raises(ValueError):step(d,u,ref,duration,100*ref.duration,True)
    assert np.array_equal(before,u)
