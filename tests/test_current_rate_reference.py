"""Independent continuum inertia operator and solver verification."""
from dataclasses import replace
import numpy as np
from scipy.special import j1,jn_zeros
from examples.current_rate_reference import CurrentRateReference
from examples.electron_reference import ElectronReference
from plasma_bench.mhd_mesh import MHDMesh


def test_variable_density_continuum_inertia_manufactured_solution_and_cell_mean():
    ref=ElectronReference();oracle=CurrentRateReference(ref,32,64)
    s,z=oracle.s[:,None],oracle.z[None];k=2*np.pi/ref.length
    exact=(1-s)**3*np.cos(k*z);ds=-3*(1-s)**2*np.cos(k*z);dss=6*(1-s)*np.cos(k*z)
    source=exact+oracle.d2/oracle.relative_density*(-4/ref.radius**2*(s*dss+2*ds-oracle.ascalar*(s*ds+exact))
        +k*k*exact-oracle.az*k*(1-s)**3*np.sin(k*z))
    value,diagnostic=oracle.solve(source)
    assert diagnostic['relative_residual']<1e-11
    assert np.linalg.norm(value-exact)<3e-11*np.linalg.norm(exact)
    mesh=MHDMesh(ref.radius,ref.length/2,13,26)
    lo=(mesh.rf[:-1]/ref.radius)**2;hi=(mesh.rf[1:]/ref.radius)**2
    expected=((1-lo)**4-(1-hi)**4)/(4*(hi-lo))
    expected=expected[:,None]*np.cos(k*mesh.z[None])*np.sinc(k*mesh.dz/(2*np.pi))
    np.testing.assert_allclose(oracle.cell_average(mesh,exact),expected,rtol=1e-10,atol=2e-15)


def test_uniform_density_bessel_inertia_eigenvalue_and_explicit_axis_probe():
    ref=replace(ElectronReference(),density_contrast=0.)
    oracle=CurrentRateReference(ref,48,64);q=jn_zeros(1,1)[0]/ref.radius;k=2*np.pi/ref.length
    r=ref.radius*np.sqrt(oracle.s[:,None])
    radial=np.divide(j1(q*r),r,out=np.full_like(r,q/2),where=r!=0)
    value=radial*np.cos(k*oracle.z[None]);exact=(1+oracle.d2*(q*q+k*k))*value;exact[-1]=value[-1]
    assert np.linalg.norm(oracle.apply(value)-exact)<3e-10*np.linalg.norm(exact)
    ref=CurrentRateReference(ElectronReference(),32,64)
    source=ref.original_source();probe=ref.original_source(axis_probe_fraction=1e-6)
    assert np.linalg.norm(probe-source)<1e-11*np.linalg.norm(source)
