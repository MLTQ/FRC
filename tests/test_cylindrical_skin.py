"""Independent continuum identities and radial induction boundary convergence."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from examples.cylindrical_skin_reference import CylindricalSkin
from examples.cylindrical_skin_measurements import RadialSkin


def test_skin_solution_satisfies_maxwell_and_diffusion():
    ref = CylindricalSkin()
    r = np.array([.13, .41, .77, .94]); h = 2e-4
    values = [ref.phasors(r+i*h) for i in (-2, -1, 0, 1, 2)]
    psi, b, current = values[2]
    first = lambda k: (values[0][k]-8*values[1][k]+8*values[3][k]-values[4][k])/(12*h)
    second = (-values[0][0]+16*values[1][0]-30*psi+16*values[3][0]-values[4][0])/(12*h*h)
    np.testing.assert_allclose(first(0)/r, b, rtol=2e-8, atol=ref.field_amplitude*1e-10)
    np.testing.assert_allclose(-first(1)/MU0, current, rtol=2e-8, atol=1e-11)
    np.testing.assert_allclose(ref.diffusivity*(second-first(0)/r), 1j*ref.omega*psi, rtol=2e-7, atol=1e-13)
    np.testing.assert_allclose(ref.diffusivity*MU0*current, -1j*ref.omega*psi/r, rtol=3e-14)


def test_skin_axis_vacuum_and_current_jump_are_explicit():
    ref = CylindricalSkin()
    r = np.array([0., 1., 1.2, 2.])
    psi, b, current = ref.phasors(r)
    assert psi[0] == 0 and current[0] == 0
    assert np.isfinite(b[0]) and abs(b[0]) < ref.field_amplitude
    np.testing.assert_allclose(b[1:], ref.field_amplitude, rtol=1e-14)
    assert current[1] != 0 and np.all(current[2:] == 0)
    np.testing.assert_allclose(psi[2:]-psi[1], .5*ref.field_amplitude*(r[2:]**2-1), rtol=2e-15)


def test_skin_cycle_power_matches_independent_volume_integral():
    ref = CylindricalSkin(frequency=7.3, diffusivity=.17, radius=.8)
    x, w = np.polynomial.legendre.leggauss(100)
    r = .5*ref.radius*(1+x)
    current = ref.phasors(r)[2]
    volume_power = np.pi*ref.diffusivity*MU0*np.sum(r*abs(current)**2*w*.5*ref.radius)
    assert ref.mean_ohmic_power_per_length() > 0
    assert abs(volume_power/ref.mean_ohmic_power_per_length()-1) < 5e-12


def test_radial_physical_current_converges_interior_and_surface():
    ref = CylindricalSkin(); errors = []
    for n in (32, 64, 128):
        s = RadialSkin(n, ref, method="physical")
        psi, _, exact = ref.fields(s.mesh.rf, .037)
        current = s.fields(psi, .037)[1]
        select = (s.mesh.rf > .2)&(s.mesh.rf < .8)
        errors.append((np.linalg.norm((current-exact)[select])/np.sqrt(select.sum()), abs(current[-1]-exact[-1])))
    errors = np.array(errors)
    assert np.all(errors[:-1, 0]/errors[1:, 0] > 3.)
    # The half-cell boundary derivative is first order pointwise, not second.
    assert np.all(errors[:-1, 1]/errors[1:, 1] > 1.6)


@pytest.mark.parametrize("method", ["physical", "standard", "strip"])
def test_evolved_skin_field_converges_to_independent_solution(method):
    errors = []
    for n in (32, 64):
        solver = RadialSkin(n, method=method)
        solver.advance(.0125)
        ref, mesh = solver.reference, solver.mesh
        exact = ref.fields(mesh.rf, solver.time)[0]
        b = solver.fields(solver.psi, solver.time)[0]
        exact_b = np.diff(exact)/(mesh.r*mesh.dr)
        errors.append(np.sqrt(np.sum((b-exact_b)**2*mesh.radial_measure)/mesh.radial_measure.sum())/ref.field_amplitude)
        assert solver.psi[0] == 0
    assert errors[1] < .012
    assert errors[0]/errors[1] > 2.
