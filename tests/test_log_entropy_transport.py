"""Continuum transport identities and independent time/space accuracy checks."""
import numpy as np
import pytest
from scipy.constants import elementary_charge as E
from scipy.linalg import expm
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.log_entropy_transport import LogEntropyTransport
from examples.log_entropy_reference import LogEntropyReference
from examples.finite_mass_manufactured import cell_average
from examples.compatible_field_norms import direct_basis


def initialize(nr=8):
    ref=LogEntropyReference();mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'])
    m=VariableCompatibleCurrent(mesh,n,ref.ion_mass)
    g=cell_average(mesh,lambda r,z:ref.fields(r,z)['g'])
    j=m.current(m.geometry.pack(np.zeros((nr+1,2*nr+1)),g))
    thermal=CompatibleLogThermal(m,gamma=ref.gamma)
    s=ref.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))['entropy']
    return ref,m,thermal,LogEntropyTransport(thermal,j),s


def test_continuum_sheared_entropy_equation_and_closed_form_integrals():
    ref=LogEntropyReference();mesh=MHDMesh(ref.radius,ref.half_length,8,16)
    r,z=mesh.centers;t=.37*ref.transit_time;h=1e-28
    f=ref.fields(r,z,t)
    st=np.imag(ref.fields(r,z,t+1j*h)['entropy'])/h
    sz=np.imag(ref.fields(r,z+1j*h,t)['entropy'])/h
    np.testing.assert_allclose(st,f['current'][...,2]*sz/(E*f['n']),rtol=3e-14,atol=1e-8)
    exact=ref.integrals()
    for time in (0.,t,ref.transit_time):
        for key in ('internal','proxy'):
            integral=np.sum(mesh.volumes*cell_average(mesh,lambda r,z:ref.fields(r,z,time)[key],order=14))
            np.testing.assert_allclose(integral,exact[key],rtol=2e-13)


def test_rk4_against_independent_physical_matrix_exponential_and_atomic_failure():
    ref,m,thermal,transport,s=initialize(6);geo=m.geometry;mesh=m.mesh
    nx,_,_=direct_basis(.5*mesh.rf**2,.5*geo.radial_points**2,5)
    nz,dz,_=direct_basis(mesh.zf,geo.axial_points,5,True)
    physical=ref.fields(geo.radial_points,np.zeros_like(geo.radial_points))
    mr=nx.T@(geo.wr[:,None]*physical['n'][:,None]*nx)
    mz=nz.T@(geo.wz[:,None]*nz)
    jr=nx.T@(geo.wr[:,None]*(physical['current'][:,2]/E)[:,None]*nx)
    az=nz.T@(geo.wz[:,None]*dz)
    operator=np.kron(np.linalg.solve(mr,jr),np.linalg.solve(mz,az))
    end=.25*ref.transit_time;target=(expm(end*operator)@s.ravel()).reshape(s.shape)
    errors=[]
    for steps in (4,8,16):
        value=s.copy()
        for _ in range(steps):value=transport.step(value,end/steps).entropy
        errors.append(np.linalg.norm(value-target)/np.linalg.norm(target-ref.background_entropy))
    assert errors[1]<.1*errors[0] and errors[2]<.1*errors[1],errors
    original=s.copy()
    with pytest.raises(ValueError,match='no floor'):transport.step(s,ref.transit_time*1e10)
    np.testing.assert_array_equal(s,original)


def test_physical_axis_transport_spatial_refinement_and_entropy_conservation():
    errors=[];axis=[]
    for nr in (8,16):
        ref,m,thermal,transport,s=initialize(nr);geo=m.geometry
        before=float(np.sum(geo.weight*m.density*thermal.scalar.sample(s)))
        end=.125*ref.transit_time;steps=2*nr
        for _ in range(steps):s=transport.step(s,end/steps).entropy
        actual=thermal.scalar.sample(s);wanted=ref.fields(geo.r,geo.z,end)['entropy']
        signal=wanted-ref.background_entropy
        error=(actual-wanted)**2;norm=signal**2
        errors.append(np.sqrt(np.sum(geo.weight*error)/np.sum(geo.weight*norm)))
        axis.append(np.sqrt(np.sum(geo.weight[:2]*error[:2])/np.sum(geo.weight[:2]*norm[:2])))
        after=float(np.sum(geo.weight*m.density*actual))
        assert abs(after-before)<2e-12*abs(before)
    assert errors[1]<.2*errors[0] and axis[1]<.2*axis[0],(errors,axis)
