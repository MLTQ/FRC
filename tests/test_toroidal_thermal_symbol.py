"""Continuum symbol cross-check through independent physical PDE differences."""
import numpy as np
from examples.electron_reference import ElectronReference
from examples.toroidal_thermal_symbol import matrices,spectrum
from examples.toroidal_symbol_difference import direct_matrix


def test_full_moving_symbol_matches_independent_physical_curl_and_species_equations():
    ref=ElectronReference();point=(.015542372881355931,-.08,1.5882496193148399)
    _,_,full,scale=matrices(ref,*point)
    expected=full*scale[None,:]/scale[:,None]
    errors=[]
    for h in (4e-6,2e-6):
        observed=direct_matrix(ref,*point,scale,spacing=h)
        errors.append(np.linalg.norm(observed-expected)/np.linalg.norm(expected))
    assert max(errors)<3e-8,errors


def test_joint_and_moving_massless_toroidal_systems_retain_the_bad_direction():
    ref=ElectronReference();s=spectrum(ref,.015542372881355931,-.08,1.5882496193148399)
    assert s['pressure_growth_per_k']>4300 and s['joint_growth_per_k']>4290
    assert max(abs(s['moving_eigenvalues'].imag))>4000


def test_uniform_density_and_pressure_gradients_recover_real_characteristics():
    ref=ElectronReference(density_contrast=0.,temperature_contrast=0.)
    for direction in np.linspace(0,2*np.pi,17):
        s=spectrum(ref,.017,.021,direction)
        assert s['pressure_growth_per_k']==0 and s['joint_growth_per_k']==0
        assert max(abs(s['moving_eigenvalues'].imag))<1e-10
