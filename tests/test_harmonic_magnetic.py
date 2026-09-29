"""Bias/current and full local energy checks for exact affine field storage."""
import numpy as np
from scipy.constants import atomic_mass,mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_current import CompatibleCurrent
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from examples.finite_mass_manufactured import cell_average


def test_uniform_bias_has_exactly_zero_current_without_discarding_field():
    mesh=MHDMesh(.035,.12,8,16);m=CompatibleCurrent(mesh,1e18,39.948*atomic_mass)
    f=HarmonicMagnetic(mesh,np.zeros((9,17)),np.zeros((8,16)),.005)
    assert np.count_nonzero(f.current(m))==0
    np.testing.assert_array_equal(f.samples(m.geometry)[...,2],.005)
    np.testing.assert_allclose(f.energy_density(m.geometry),.005**2/(2*MU0),rtol=2e-15)
    assert f.energy_components(m)['increment']==0.


def test_weak_polynomial_current_survives_bias_and_local_cross_energy():
    mesh=MHDMesh(.7,.6,8,12);m=CompatibleCurrent(mesh,1e18,39.948*atomic_mass);geo=m.geometry
    x=.5*mesh.corners[0]**2;X=.5*mesh.radius**2
    psi=x*(1-x/X)**3;psi[[0,-1]]=0.
    g=cell_average(mesh,lambda r,z:(1-.5*r*r/X)**2)
    y=.5*geo.r**2/X
    j=np.stack((np.zeros_like(y),geo.r*(6-18*y+12*y*y)/(X*MU0),2*(1-4*y+3*y*y)/MU0),axis=-1)
    perturb=np.stack((np.zeros_like(y),geo.r*(1-y)**2,1-6*y+9*y*y-4*y**3),axis=-1)
    for amplitude in (1.,1e-4,1e-8,1e-12):
        field=HarmonicMagnetic(mesh,amplitude*psi,amplitude*g,.4)
        actual=geo.sample_electric(field.current(m))/amplitude
        np.testing.assert_allclose(actual,j,rtol=3e-9,atol=3e-6)
        energies=field.energy_components(m)
        expected=np.sum(geo.weight*np.sum(perturb*perturb,axis=-1))/(2*MU0)
        np.testing.assert_allclose(energies['increment']/amplitude**2,expected,rtol=2e-13)
        full=amplitude*perturb;full[...,2]+=.4
        cells=np.sum(geo.weight*np.sum(full*full,axis=-1),axis=(2,3))/(2*MU0*mesh.volumes)
        np.testing.assert_allclose(field.energy_density(geo),cells,rtol=3e-13)
        np.testing.assert_allclose(np.sum(cells*mesh.volumes),sum(energies.values()),rtol=3e-13)
