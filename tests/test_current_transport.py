import numpy as np
import pytest
from examples.regular_coupling_reference import initialize
from examples.current_transport_study import run_case
from plasma_bench.current_transport import CurrentTransport


def test_actual_polynomial_current_and_covariant_velocity_derivatives():
    ref,d,u,m,field=initialize(dict(elements=4));g=m.geometry;t=CurrentTransport(m);y=(g.r/ref.radius)**2
    values=np.stack((g.r*(1+.2*y),g.r*(1-y),1+.3*y),axis=-1);coeff=g.project_electric(values);grad=t.current.gradient(coeff)
    wanted=np.zeros_like(grad);wanted[...,0,0]=1+.6*y;wanted[...,1,0]=1-3*y;wanted[...,2,0]=.6*g.r/ref.radius**2
    wanted[...,0,1]=-(1-y);wanted[...,1,1]=1+.2*y
    assert np.linalg.norm(grad-wanted)/np.linalg.norm(wanted)<3e-12
    n=u[...,0];u[...,1]=.12*n;u[...,2]=.2*n;u[0,0,:,1:3]=0
    k=np.pi/ref.half_length;u[...,3]=n*(3*np.cos(k*d.z)+2)
    vg=t.fluid.velocity_gradient(u);exact=np.zeros_like(vg);exact[...,0,0]=exact[...,1,1]=.12;exact[...,0,1]=-.2;exact[...,1,0]=.2;exact[...,2,2]=-3*k*np.sin(k*g.z)
    assert np.linalg.norm(vg-exact)/np.linalg.norm(exact)<3e-12
    stress=np.zeros_like(vg);stress[...,0,0]=1+.1*y;stress[...,1,1]=.7;stress[...,2,2]=2+.1*y
    stress[...,0,2]=stress[...,2,0]=g.r*np.sin(k*g.z)
    rate=t.fluid.stress_rate(u,stress);work=np.sum(d.volume[...,None]*d.fields(u)['energy_gradient']*rate)
    wanted=np.sum(g.weight*np.sum(vg*stress,axis=(-2,-1)))
    assert abs(work-wanted)<1e-12*max(abs(work),abs(wanted))


def test_original_local_terms_refine_and_face_omission_is_inaccurate():
    coarse,_=run_case(dict(elements=4));fine,_=run_case(dict(elements=8))
    for key in ('advection','stretching','self_drift'):
        assert coarse['electric_relative_L2'][key]/fine['electric_relative_L2'][key]>5
    assert np.all(np.array(coarse['force_relative_L2'])/fine['force_relative_L2']>10)
    assert fine['advection_without_faces']>fine['electric_relative_L2']['advection']
    assert np.all(np.array(fine['force_without_faces'])>10*np.array(fine['force_relative_L2']))


def test_actual_current_density_exchange_and_drift_stress_close_without_heat():
    r,_=run_case(dict(elements=8))
    for key in ('advection_work_backward','stress_work_backward','self_drift_work_backward','axial_force_scaled','angular_torque_scaled'):
        assert r[key]<2e-12
    assert abs(r['density_exchange'])>1e-4 and abs(r['stress_work'])>1e-4
    assert r['total_energy_backward']<1e-10 and r['residual']<1e-10


def test_transport_rejects_mismatched_or_invalid_state_atomically():
    ref,d,u,m,field=initialize(dict(elements=2));t=CurrentTransport(m);current=field.current(m);before=u.copy();nt=d.derivative(u)[0][...,0]
    bad=u.copy();bad[...,0]*=1.01
    with pytest.raises(ValueError):t.pair(bad,current,nt)
    bad=u.copy();bad[0,0,:,1]=1.
    with pytest.raises(ValueError):t.pair(bad,current,nt)
    assert np.array_equal(u,before)
