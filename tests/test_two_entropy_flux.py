import numpy as np
import pytest
from scipy.constants import elementary_charge as E
from plasma_bench.two_entropy_flux import TwoEntropyFlux
from plasma_bench.two_entropy_axial import TwoEntropyAxial
from examples.log_fluid_reference import LogFluidReference


def initial(solver,ref,time=0.):
    f=ref.fields(np.zeros_like(solver.z),solver.z,time)
    return np.concatenate((f['n'][:,None],f['n'][:,None]*f['velocity'],f['n'][:,None]*np.stack((f['electron'],f['ion']),axis=-1)),axis=-1)


def energy(state,ref):
    n=state[...,0];s=state[...,4:]/n[...,None];v=state[...,1:4]/n[...,None]
    h=E*np.exp(s)*(n[...,None]/1e20)**(np.array([ref.gamma_e,ref.gamma_i])-1)/(np.array([ref.gamma_e,ref.gamma_i])-1)
    return .5*ref.mass*n*np.sum(v*v,axis=-1)+n*np.sum(h,axis=-1)


def test_physical_gradient_and_pair_flux_identity_across_independent_states():
    ref=LogFluidReference();model=TwoEntropyFlux(ref.mass);rng=np.random.default_rng(741)
    n=ref.number_density*np.exp(rng.uniform(-2,2,(80,2)));v=ref.sound*rng.normal(size=(80,2,3))
    s=rng.uniform(-1,3,(80,2,2));u=np.concatenate((n[...,None],n[...,None]*v,n[...,None]*s),axis=-1)
    f=model.fields(u);direction=u*rng.normal(size=u.shape);oracle=np.imag(energy(u+1j*1e-28*direction,ref))/1e-28
    assert np.max(abs(np.sum(f['energy_gradient']*direction,axis=-1)-oracle)/np.maximum(abs(oracle),1e-25))<5e-12
    for axis in (0,1,2):
        flux,g=model.pair(u[:,0],u[:,1],axis);reverse,rg=model.pair(u[:,1],u[:,0],axis)
        assert np.array_equal(flux,reverse) and np.array_equal(g,rg)
        w=f['energy_gradient'];psi=np.sum(f['pressure'],axis=-1)*f['velocity'][...,axis]
        left=np.sum((w[:,1]-w[:,0])*flux,axis=-1);right=psi[:,1]-psi[:,0]
        scale=np.sum(abs((w[:,1]-w[:,0])*flux),axis=-1)+abs(right)
        assert np.max(abs(left-right)/scale)<3e-13
        same,same_g=model.pair(u[:,0],u[:,0],axis)
        assert np.allclose(same[:,0],u[:,0,1+axis],rtol=2e-15)
        assert np.allclose(same_g,(f['energy'][:,0]+np.sum(f['pressure'][:,0],axis=-1))*f['velocity'][:,0,axis],rtol=3e-14)


def test_original_contact_rates_refine_and_spatial_energy_is_conservative():
    errors=[]
    for points in (16,32,64):
        ref=LogFluidReference();solver=TwoEntropyAxial(points,2*ref.half_length,ref.mass)
        u=initial(solver,ref);rate,d=solver.derivative(u);f=ref.fields(np.zeros(points),solver.z)
        _,_,exact=ref.derivatives(np.zeros(points),solver.z)
        pressure=[]
        for a,key,gamma in ((0,'pe',ref.gamma_e),(1,'pi',ref.gamma_i)):
            s=u[:,4+a]/u[:,0];st=(rate[:,4+a]-s*rate[:,0])/u[:,0]
            pressure.append(f[key]*(st+gamma*rate[:,0]/u[:,0]))
        errors.append(np.linalg.norm(pressure[1]-exact['pi'])/np.linalg.norm(exact['pi']))
        assert np.max(abs(np.sum(rate,axis=0))/np.maximum(np.sum(abs(rate),axis=0),1e-30))<2e-12
        assert abs(np.sum(d['physical_energy_rate']))/np.sum(abs(d['physical_energy_rate']))<5e-12
        assert np.linalg.norm(d['physical_energy_rate']-d['energy_flux_rate'])/np.linalg.norm(d['physical_energy_rate'])<5e-12
        # Uniform velocity must not acquire a force from a global energy fix.
        vt=(rate[:,1:4]-f['velocity']*rate[:,0,None])/f['n'][:,None]
        assert np.max(abs(vt))/ref.sound**2*ref.half_length<3e-13
    assert all(a/b>30 for a,b in zip(errors[:-1],errors[1:]))


def test_nonlinear_wave_accuracy_and_atomic_rejection():
    errors=[]
    for points in (16,32):
        ref=LogFluidReference(kind='acoustic',gamma_i=5/3);solver=TwoEntropyAxial(points,2*ref.half_length,ref.mass,gammas=(ref.gamma_e,ref.gamma_i))
        u=initial(solver,ref);duration=.2*ref.duration
        for _ in range(40):u=solver.step(u,duration/40)
        exact=ref.fields(np.zeros(points),solver.z,duration)['n'];errors.append(np.linalg.norm(u[:,0]-exact)/np.linalg.norm(exact-ref.number_density))
    assert errors[0]/errors[1]>30
    before=u.copy()
    with pytest.raises(ValueError):solver.step(u,100*ref.duration)
    assert np.array_equal(before,u)
