"""Independent Maxwell/Ohm, interface and whole-space checks of the decay target."""
import numpy as np
from scipy.constants import mu_0 as MU0
from examples.spherical_decay_reference import SphericalDecay


def test_spherical_decay_matches_maxwell_faraday_and_ohm():
    ref=SphericalDecay();r=np.array([.12,.43,.71,1.4]);z=np.array([.3,-.17,.2,.2])
    t=.13;h=1e-5;f=ref.fields(r,z,t)
    dr=(ref.fields(r+h,z,t)["psi"]-ref.fields(r-h,z,t)["psi"])/(2*h*r)
    dz=-(ref.fields(r,z+h,t)["psi"]-ref.fields(r,z-h,t)["psi"])/(2*h*r)
    current=((ref.fields(r,z+h,t)["Br"]-ref.fields(r,z-h,t)["Br"])-
             (ref.fields(r+h,z,t)["Bz"]-ref.fields(r-h,z,t)["Bz"]))/ (2*h*MU0)
    electric=-(ref.fields(r,z,t+h)["psi"]-ref.fields(r,z,t-h)["psi"])/(2*h*r)
    np.testing.assert_allclose(dr,f["Bz"],atol=2e-12,rtol=2e-8)
    np.testing.assert_allclose(dz,f["Br"],atol=2e-12,rtol=2e-8)
    np.testing.assert_allclose(current,f["Jtheta"],atol=2e-6,rtol=2e-8)
    np.testing.assert_allclose(electric,f["Etheta"],atol=1e-12,rtol=2e-8)
    np.testing.assert_allclose(f["Etheta"][:3],MU0*ref.diffusivity*f["Jtheta"][:3],rtol=2e-14)
    assert f["Jtheta"][-1]==0 and f["Etheta"][-1]!=0


def test_spherical_decay_axis_interface_and_analytic_whole_space_energy():
    nodes,weights=np.polynomial.legendre.leggauss(72);x=(nodes+1)/2;w=weights/2
    for mode in (1,2):
        ref=SphericalDecay(mode=mode);t=.07
        a,ap=ref.radial(np.array([0.,1.]),t)
        assert ap[0]==0
        np.testing.assert_allclose(ap[1],-3*a[1],atol=1e-17)
        assert abs(ref.fields(0.,0.,t)["Bz"]-2*ref.central_a*np.exp(-ref.diffusivity*ref.wavenumber**2*t)) < 1e-17
        energy=0.;heat=0.
        for s,radial_weight in ((x,w*x*x),(1/x,w/x**4)):
            r=s[:,None]*np.sqrt(1-nodes[None,:]**2);z=s[:,None]*nodes[None,:]
            f=ref.fields(r,z,t)
            energy+=np.pi/MU0*np.sum(radial_weight[:,None]*weights[None,:]*(f["Br"]**2+f["Bz"]**2))
            heat+=2*np.pi*np.sum(radial_weight[:,None]*weights[None,:]*ref.heat_density(r,z,t))
        exact=ref.integrals(t)
        assert abs(energy/exact["induced_energy_j"]-1) < 2e-13
        assert abs(heat/exact["ohmic_heat_j"]-1) < 2e-13
