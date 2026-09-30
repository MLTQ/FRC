import numpy as np
import pytest
from plasma_bench.gauss_sbp import GaussSBP
from plasma_bench.cylindrical_entropy_fluid import CylindricalEntropyFluid
from plasma_bench.cylindrical_entropy_state import CylindricalEntropyState
from plasma_bench.thermal_means import logarithmic_mean
from examples.log_fluid_reference import LogFluidReference
from examples.radial_work_reference import RadialWorkReference
from examples.cylindrical_entropy_reference import initial,physical,rates,regular_sample


def solver(ref,elements=4,points=32,degree=5):
    return CylindricalEntropyFluid(ref.radius,2*ref.half_length,elements,points,ref.mass,degree=degree,gammas=(ref.gamma_e,ref.gamma_i))


def test_gauss_polynomial_integrals_derivatives_and_boundary_traces():
    for degree in (1,3,5,7):
        b=GaussSBP(degree)
        for p in range(degree+1):
            exact=p*b.nodes**(p-1) if p else np.zeros_like(b.nodes)
            assert np.max(abs(b.derivative@b.nodes**p-exact))<2e-12
            assert np.max(abs(b.trace@b.nodes**p-np.array([(-1.)**p,1.])))<2e-13
        for p in range(2*degree+2):
            assert abs(b.weights@b.nodes**p-(2/(p+1) if p%2==0 else 0))<2e-14


def test_actual_cylindrical_energy_gradient_inverse_and_geometric_pair_jump():
    ref=RadialWorkReference();d=solver(ref);u=initial(d,ref);f=d.fields(u)
    rng=np.random.default_rng(739);direction=u*rng.normal(size=u.shape)
    wanted=np.imag(physical(u+1j*1e-28*direction,d.r,ref)['energy'])/1e-28
    assert np.linalg.norm(np.sum(f['energy_gradient']*direction,axis=-1)-wanted)/np.linalg.norm(wanted)<3e-14
    back,_=d.state_map.from_energy_variables(f['energy_gradient'],d.r)
    assert np.max(abs(back-f['physical_state'])/np.maximum(abs(f['physical_state']),1))<2e-13
    trace=np.einsum('bq,eqzk->ebzk',d.sbp.trace,f['energy_gradient'])
    boundary,bf=d.state_map.from_energy_variables(trace,d.boundary_radius)
    assert np.array_equal(boundary[0,0,:,2],np.zeros(d.nz))
    # Different radii distinguish actual angular momentum from Cartesian mt.
    left=f['physical_state'][:,0];right=f['physical_state'][:,-1];rl=d.r[:,0];rr=d.r[:,-1]
    flux=d.state_map.pair(left,right,rl,rr,0)
    jump=f['energy_gradient'][:,-1]-f['energy_gradient'][:,0]
    psi=np.sum(f['pressure'],axis=-1)*f['velocity'][...,0]
    gamma=logarithmic_mean(left[...,0],right[...,0])*.5*(f['velocity'][:,0,:,0]+f['velocity'][:,-1,:,0])
    extra=ref.mass*gamma*(rr-rl)*.5*(f['velocity'][:,0,:,1]+f['velocity'][:,-1,:,1])*.5*(f['omega'][:,0]+f['omega'][:,-1])
    lhs=np.sum(jump*flux,axis=-1);rhs=psi[:,-1]-psi[:,0]-extra
    assert np.max(abs(lhs-rhs))/np.max(abs(rhs))<3e-13


def test_original_2d_axis_rates_converge_and_element_energy_is_conservative():
    ref=RadialWorkReference();errors=[]
    for ne in (4,8,16):
        d=solver(ref,ne,64);u=initial(d,ref);du,_=d.derivative(u);actual=rates(u,du,d.r,ref);wanted=ref.rates(d.r,d.z)
        error=actual['velocity']-wanted['velocity']
        errors.append(np.sqrt(np.sum(d.volume[0]*error[0,...,0]**2)/np.sum(d.volume[0]*wanted['velocity'][0,...,0]**2)))
        summed=np.sum(d.volume[...,None]*du,axis=(0,1,2));scale=np.sum(d.volume[...,None]*abs(du),axis=(0,1,2))
        assert np.max(abs(summed[[0,2,3,4,5]])/scale[[0,2,3,4,5]])<1e-11
        assert abs(np.sum(d.volume*actual['energy']))/np.sum(d.volume*abs(actual['energy']))<2e-12
        # Independent physical-energy chain rule balances the shared boundary
        # energy flux for each individual radial element, not only globally.
        _,b,_,_,_=d.extended(u);r=d.boundary_radius
        fleft=d.state_map.thermal.fields(b[:-1,1]);fright=d.state_map.thermal.fields(b[1:,0])
        flux,_=d.state_map.thermal.pair(b[:-1,1],b[1:,0],0)
        psi_l=np.sum(fleft['pressure'],axis=-1)*fleft['velocity'][...,0]
        psi_r=np.sum(fright['pressure'],axis=-1)*fright['velocity'][...,0]
        g=.5*np.sum((fleft['energy_gradient']+fright['energy_gradient'])*flux,axis=-1)-.5*(psi_l+psi_r)
        interface=2*np.pi*d.dz*np.sum(r[:-1,1]*g,axis=1);outward=np.diff(np.r_[0.,interface,0.])
        integrated=np.sum(d.volume*actual['energy'],axis=(1,2))
        assert np.max(abs(integrated+outward))/np.max(abs(integrated))<2e-11
    assert errors[0]/errors[1]>5 and errors[1]/errors[2]>10


def test_regular_physical_reconstruction_has_exact_axis_parity_and_polynomial_energy():
    ref=LogFluidReference();d=solver(ref,3,16);r,z=np.broadcast_arrays(d.r,d.z)
    n=ref.number_density*(1+.2*(r/ref.radius)**2);u=np.zeros(d.shape);u[...,0]=n
    u[...,1]=n*.1*ref.sound*r/ref.radius;u[...,2]=n*.2*ref.sound*r*r/ref.radius
    u[...,3]=n*.3*ref.sound;u[...,4]=n*.2;u[...,5]=n*.4
    sample,rr,zz,vol=regular_sample(d,u,24,2);nn=ref.number_density*(1+.2*(rr/ref.radius)**2)
    assert np.max(abs(sample[...,0]-nn)/nn)<2e-13
    v=physical(sample,rr,ref)['velocity']
    assert np.max(abs(v[...,0]-.1*ref.sound*rr/ref.radius))/ref.sound<2e-13
    assert np.max(abs(v[...,1]-.2*ref.sound*rr/ref.radius))/ref.sound<2e-13
    # Integrating a constant recovers the exact cylindrical volume.
    assert abs(np.sum(vol)/(np.pi*ref.radius**2*d.length)-1)<2e-14


def test_uniform_translation_has_no_spurious_radial_force_and_rejection_is_atomic():
    ref=LogFluidReference();d=solver(ref,3,32);u=initial(d,ref);du,_=d.derivative(u)
    actual=rates(u,du,d.r,ref)
    assert np.max(abs(actual['velocity']))*ref.radius/ref.sound**2<2e-11
    before=u.copy()
    with pytest.raises(ValueError):d.step(u,100*ref.duration)
    assert np.array_equal(u,before)


def test_radial_acoustic_spectrum_matches_bessel_roots_and_original_nonlinear_rhs():
    from examples.cylindrical_entropy_spectrum import calculate
    errors=[]
    for ne in (2,4):
        record,_=calculate(dict(elements=ne,degree=3))
        errors.append(abs(record['first_three_frequency_relative_errors'][0]))
        assert record['energy_adjoint_defect']<1e-12
        assert record['maximum_real_scaled']<1e-11
        assert record['nonlinear_difference_controls'][1]['relative_error']<2e-9
        temporal=record['rk4_controls']
        assert temporal[0]['relative_error']/temporal[1]['relative_error']>15
    assert errors[0]/errors[1]>40
