import numpy as np
import pytest
from scipy.constants import elementary_charge as E
from plasma_bench.regular_field_map import periodic,RegularDensity
from plasma_bench.regular_current import RegularCurrent
from examples.regular_coupling_reference import initialize
from examples.regular_coupling_study import run_case
from examples.regular_coupling_norms import measure


def test_periodic_actual_interpolant_including_nyquist_and_regular_density():
    for n in (15,16):
        x=np.arange(n)/n;target=np.linspace(-.2,1.3,71);p=periodic(x,target,1.)
        for k in range(n//2+1):assert np.max(abs(p@np.cos(2*np.pi*k*x)-np.cos(2*np.pi*k*target)))<2e-14
    ref,d,u,m,_=initialize(dict(elements=2));g=m.geometry
    assert np.max(abs(m.density/ref.fields(g.r,g.z)['n']-1))<2e-14
    before=u.copy()
    with pytest.raises(ValueError):RegularDensity(m.mapping,-u[...,0])
    with pytest.raises(ValueError):RegularCurrent(m.mapping,u[...,0],2*ref.ion_mass)
    assert np.array_equal(before,u)


def test_independent_physical_forces_refine_and_pair_energy_momentum_torque():
    errors=[];axis=[]
    for ne in (2,4,8):
        r,_=run_case(dict(elements=ne));errors.append(np.array(r['force_relative_L2']));axis.append(r['axis_force_relative_L2'][2])
        for key in ('force_work_backward','axial_force_relative','angular_torque_relative','density_gradient_relative','ideal_work_backward','hall_work_backward'):
            assert r[key]<2e-12
        assert r['induction_energy_backward']<1e-10 and r['exchange_balance_backward']<1e-10
    assert np.all(errors[0]/errors[1]>15) and np.all(errors[1]/errors[2]>15)
    assert axis[0]/axis[1]>100 and axis[1]/axis[2]>100


def test_original_induction_and_lorentz_rates_refine_without_concealing_exchange():
    coarse,_=run_case(dict(elements=4));fine,_=run_case(dict(elements=8))
    for key in ('current_error','magnetic_rate_error','current_rate_error'):
        assert coarse[key]/fine[key]>20
    assert np.all(np.array(coarse['lorentz_force_error'])/fine['lorentz_force_error']>10)
    assert abs(fine['unresolved_density_exchange'])>1e-4
    assert fine['residual']<=1e-10 and fine['coupled_residual']<=1e-10


def test_independent_root_fft_energy_audit_and_atomic_pair_rejection():
    spec=dict(elements=4);r,a=run_case(spec);ref,d,u,m,field=initialize(spec)
    measured=measure(d,m.mesh,ref,a,order=16);j=m.geometry.sample_electric(a['current'])
    energy=float(np.sum(m.geometry.weight*.5*m.reduced_mass*np.sum(j*j,axis=-1)/(E**2*m.density)))
    assert abs(measured['current_energy']/energy-1)<3e-12
    assert abs(measured['density_exchange']/r['unresolved_density_exchange']-1)<3e-10
    bulk=float(np.sum(d.volume[...,None]*d.fields(u)['energy_gradient']*a['lorentz_rate']))
    assert abs(measured['lorentz_work']/bulk-1)<3e-12
    before=u.copy();wrong=u.copy();wrong[...,0]*=1.01
    with pytest.raises(ValueError):m.electromagnetic_pair(wrong,field,a['fluid_rate'][...,0])
    assert np.array_equal(u,before)
