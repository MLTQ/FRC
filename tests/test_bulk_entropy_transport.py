"""Independent changing-density characteristics and atomic RK4 verification."""
import numpy as np
import pytest
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.bulk_entropy_transport import BulkEntropyState,PrescribedBulkTransport
from examples.bulk_compression_reference import BulkCompressionReference
from examples.finite_mass_manufactured import cell_average
from examples.physical_scalar_sampler import PhysicalScalarSampler


def initial(ref,nr=6):
    mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr);geometry=CompatibleGeometry(mesh,quadrature_order=12)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=14)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=ref.fields(r,z)
    state=BulkEntropyState(n,f['electron'],f['ion'])
    return PrescribedBulkTransport(geometry,ref.ion_mass,ref.momentum,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i),state


def test_reference_parent_continuity_and_both_adiabatic_pressure_equations():
    for kind in ('compression','translation'):
        ref=BulkCompressionReference(kind=kind);r,z=np.meshgrid(np.linspace(0,ref.radius,5),np.linspace(-.113,.115,9),indexing='ij')
        t=.47*ref.duration;h=1e-28;f=ref.fields(r,z,t)
        ft=ref.fields(r,z,t+1j*h);fz=ref.fields(r,z+1j*h,t)
        dv=np.imag(ref.velocity(z+1j*h))/h
        np.testing.assert_allclose(np.imag(ft['n'])/h,-f['velocity']*np.imag(fz['n'])/h-f['n']*dv,rtol=2e-13,atol=5e7)
        for label,gamma in (('electron',ref.gamma_e),('ion',ref.gamma_i)):
            pt=np.imag(ft[label+'_pressure'])/h
            rhs=-f['velocity']*np.imag(fz[label+'_pressure'])/h-gamma*f[label+'_pressure']*dv
            np.testing.assert_allclose(pt,rhs,rtol=2e-13,atol=1e-7)


def test_changing_density_time_order_work_and_atomic_rejection():
    ref=BulkCompressionReference();solver,state=initial(ref);end=.3*ref.duration
    results=[];work_errors=[]
    def internal(value):
        model,_,fields=solver.evaluate(value)
        return np.array([np.sum(model.geometry.weight*f['internal']) for f in fields])
    initial_energy=internal(state)
    for steps in (2,4,8):
        value=state;work=np.zeros(2)
        for step in range(steps):
            result=solver.step(value,step*end/steps,end/steps);value=result.state;work+=result.thermal_work
        results.append(value);work_errors.append(np.linalg.norm(internal(value)-initial_energy-work))
    d1=np.linalg.norm(results[0].electron-results[1].electron);d2=np.linalg.norm(results[1].electron-results[2].electron)
    assert 12<d1/d2<20,(d1,d2)
    assert work_errors[1]<.12*work_errors[0] and work_errors[2]<.12*work_errors[1],work_errors
    copy=[getattr(state,k).copy() for k in ('density','electron','ion')]
    with pytest.raises(ValueError):solver.step(state,0,1e9*ref.duration)
    for k,value in zip(('density','electron','ion'),copy):np.testing.assert_array_equal(getattr(state,k),value)
