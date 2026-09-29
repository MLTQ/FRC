"""Independent derivative checks for the non-Beltrami nonlinear inertia fixture."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import electron_mass as ME,atomic_mass,elementary_charge as E
from examples.electron_reference import ElectronReference
from examples.current_inertia_reference import fields,g_rate


@pytest.mark.parametrize('contrast',(0.,.4))
def test_full_convective_acceleration_gradient_and_nonzero_curl(contrast):
    ref=replace(ElectronReference(),density_contrast=contrast)
    mi=39.948*atomic_mass;mu=ME*mi/(ME+mi);delta=(mi-ME)/(mi+ME)
    r,z=np.meshgrid(np.linspace(.001,.034,13),np.linspace(-.11,.11,17),indexing='ij')
    value=fields(ref,r,z,mu,delta);original=ref.fields(r,z)
    np.testing.assert_allclose(value['drift'][...,0]*original['density']*E,original['jr'],rtol=5e-15,atol=1e-15)
    np.testing.assert_allclose(value['drift'][...,2]*original['density']*E,original['jz'],rtol=5e-15,atol=1e-15)
    np.testing.assert_allclose(value['full_convective_electric'],value['reduced_convective_electric']-delta/E*value['kinetic_gradient'],rtol=2e-12,atol=1e-18)
    rate=g_rate(ref,r,z,mu,delta);full=g_rate(ref,r,z,mu,delta,full=True)
    assert np.linalg.norm(rate)>1.
    assert np.linalg.norm(rate-full)<2e-14*np.linalg.norm(rate)
    h=1e-6
    def derivative(axis):
        v=[fields(ref,r+i*h if axis==0 else r,z+i*h if axis==1 else z,mu,delta)['full_convective_electric'] for i in (-2,-1,1,2)]
        return (v[0]-8*v[1]+8*v[2]-v[3])/(12*h)
    finite=(derivative(0)[...,2]-derivative(1)[...,0])/r
    assert np.linalg.norm(finite-rate)<2e-8*np.linalg.norm(rate)
