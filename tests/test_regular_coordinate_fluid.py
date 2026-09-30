import numpy as np
import pytest
from plasma_bench.lobatto_sbp import LobattoSBP
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from examples.regular_euler_reference import RegularEulerReference,initial,energy
from examples.regular_euler_sampling import basis
from examples.regular_euler_spectrum import calculate


def test_lobatto_integrals_and_positive_axis_velocity_physical_metric():
    for p in (3,5,7):
        b=LobattoSBP(p)
        for k in range(2*p):assert abs(b.weights@b.nodes**k-(2/(k+1) if k%2==0 else 0))<5e-14
        x=.5*(b.nodes+1);w=.5*b.weights;q,wq=np.polynomial.legendre.leggauss(2*p+3);r=.5*(q+1)
        v=basis(x[1:],r);gram=v.T@((.5*wq*2*r)[:,None]*v)
        assert np.linalg.norm(gram-np.diag((2*x*w)[1:]))/np.linalg.norm(gram)<5e-13
        assert np.min(np.linalg.eigvalsh(gram))>0


def test_original_regular_parent_rates_and_axis_refine_with_true_energy_and_angular_balance():
    ref=RegularEulerReference();errors=[];axis=[]
    for ne in (2,4,8):
        d=RegularCoordinateFluid(ref.radius,2*ref.half_length,ne,128,ref.mass,degree=3);u=initial(d,ref);du,diag=d.derivative(u);wanted=ref.parent_rate(d.x,d.z)
        wanted[0,0,:,1:3]=0;diff=du-wanted
        errors.append(np.sqrt(np.sum(d.volume*diff[...,4]**2)/np.sum(d.volume*wanted[...,4]**2)))
        axis.append(np.sqrt(np.sum(d.volume[0,1:]*diff[0,1:,:,1]**2)/np.sum(d.volume[0,1:]*wanted[0,1:,:,1]**2)))
        independent=np.imag(energy(u+1j*1e-28*du,d.x,ref))/1e-28
        assert np.linalg.norm(independent-diag['physical_energy_rate'])/np.linalg.norm(independent)<3e-14
        assert abs(np.sum(d.volume*independent))/np.sum(d.volume*abs(independent))<2e-12
        angular=2*d.x*du[...,2]
        assert abs(np.sum(d.volume*angular))/np.sum(d.volume*abs(angular))<2e-12
    assert all(a/b>5 for a,b in zip(errors[:-1],errors[1:]))
    assert all(a/b>40 for a,b in zip(axis[:-1],axis[1:]))


def test_analytic_bessel_frequency_converges_with_positive_acoustic_metric():
    errors=[]
    for n in (2,4,8):
        record,a=calculate(n,3);errors.append(abs(record['frequency_errors'][0]))
        assert record['maximum_real']<2e-11
        assert record['energy_adjoint_defect']<5e-12
        assert np.min(np.diag(a['metric']))>0
        assert record['nonlinear_controls'][1]['relative_error']<1e-8
        assert record['rk4_controls'][0]['relative_error']/record['rk4_controls'][1]['relative_error']>14
    assert errors[0]/errors[1]>20 and errors[1]/errors[2]>20


def test_inactive_axis_dofs_and_rejected_step_do_not_repair_state():
    ref=RegularEulerReference();d=RegularCoordinateFluid(ref.radius,2*ref.half_length,2,16,ref.mass);u=initial(d,ref);before=u.copy()
    with pytest.raises(ValueError):d.step(u,100*ref.duration)
    assert np.array_equal(before,u)
    invalid=u.copy();invalid[0,0,:,1]=1.
    with pytest.raises(ValueError):d.fields(invalid)
