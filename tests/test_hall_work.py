"""Algebraic Hall work and independent spatial consistency gates."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.toroidal_state import radius_squared
from plasma_bench.hall_curl import HallCurl
from plasma_bench.hall_induction import hall_induction
from plasma_bench.hall_poynting import local_energy_audit
from examples.hall_manufactured import HallManufactured
from examples.hall_work_snapshots import snapshot


def magnetic_energy(mesh, psi, g):
    br, bz = mesh.cell_fields(psi)
    return np.sum(mesh.volumes*(br*br+bz*bz+radius_squared(mesh)*g*g))/(2*MU0)


def test_compatible_current_differentiates_the_actual_magnetic_energy():
    mesh = MHDMesh(.7, .8, 9, 12)
    geometry = HallCurl(mesh)
    rng = np.random.default_rng(2419)
    psi = rng.normal(size=(mesh.nr+1, mesh.nz+1))*1e-4
    psi[[0, -1]] = 0.; psi[:, -1] = psi[:, 0]
    g = rng.normal(size=(mesh.nr, mesh.nz))*1e-3
    electric = rng.normal(size=geometry.offsets[-1])
    dp, dg = geometry.unpack(geometry.faraday @ electric)
    exact = np.dot(geometry.magnetic_gradient(psi, g), geometry.faraday @ electric)
    # A quadratic has no central-difference truncation error. Use a resolvable
    # perturbation rather than subtracting nearly identical SI energies.
    eps = .1*np.linalg.norm(geometry.pack(psi, g))/np.linalg.norm(geometry.pack(dp, dg))
    measured = (magnetic_energy(mesh, psi+eps*dp, g+eps*dg)-magnetic_energy(mesh, psi-eps*dp, g-eps*dg))/(2*eps)
    assert measured == pytest.approx(exact, rel=2e-8)
    work = -np.sum(geometry.weight*np.sum(geometry.current(psi, g)*geometry.sample_electric(electric), axis=-1))
    assert work == pytest.approx(exact, rel=2e-14)


@pytest.mark.parametrize('density_kind', ['scalar', 'cell', 'quadrature'])
def test_nonlinear_hall_cancels_work_and_preserves_ct_without_heat(density_kind):
    reference = HallManufactured()
    mesh = MHDMesh(reference.radius, reference.length/2, 11, 18)
    geometry = HallCurl(mesh)
    initial = snapshot(reference, mesh)
    density = {'scalar': reference.number_density, 'cell': initial['density'],
               'quadrature': reference.fields(geometry.r, geometry.z)['number_density']}[density_kind]
    result = hall_induction(geometry, initial['psi'], initial['g'], density)
    assert abs(result.magnetic_work_w)/result.absolute_work_scale_w < 2e-15
    assert np.max(abs(mesh.divergence(result.psi_rate)))*mesh.dr/(np.max(abs(geometry.magnetic_samples(result.psi_rate, result.g_rate)))) < 2e-13
    assert_allclose(result.psi_rate[[0, -1]], 0., atol=0.)
    assert_allclose(result.psi_rate[:, 0], result.psi_rate[:, -1], atol=0.)
    magnetic, total, sr, sz = local_energy_audit(geometry, initial['psi'], initial['g'], result)
    assert_allclose(sr[[0, -1]], 0., atol=0.)
    assert_allclose(sz[:, 0], sz[:, -1], atol=0.)
    scale = np.sum(mesh.volumes*abs(magnetic))
    assert abs(np.sum(mesh.volumes*magnetic))/scale < 1e-14
    assert abs(np.sum(mesh.volumes*total))/scale < 1e-14
    # Do not turn a local mismatch into a claimed Joule heat or an energy repair.
    assert np.linalg.norm(total-magnetic) > 1e-4*np.linalg.norm(magnetic)


def test_uniform_field_has_no_hall_evolution():
    mesh = MHDMesh(1., 1., 8, 8)
    geometry = HallCurl(mesh)
    r, _ = mesh.corners
    result = hall_induction(geometry, .5*r*r, np.zeros((8, 8)), 1e20)
    assert_allclose(result.psi_rate, 0., atol=2e-8)
    assert_allclose(result.g_rate, 0., atol=2e-8)


def test_independent_nonlinear_induction_refines_but_is_not_certified_at_axis():
    reference = HallManufactured()
    errors, local = [], []
    for n in (8, 16, 32):
        mesh = MHDMesh(reference.radius, reference.length/2, n, 2*n)
        geometry = HallCurl(mesh)
        initial = snapshot(reference, mesh)
        exact = reference.fields(geometry.r, geometry.z)
        result = hall_induction(geometry, initial['psi'], initial['g'], exact['number_density'])
        actual = geometry.magnetic_samples(result.psi_rate, result.g_rate)
        errors.append(np.sqrt(np.sum(geometry.weight[..., None]*(actual-exact['induction'])**2)/np.sum(geometry.weight[..., None]*exact['induction']**2)))
        magnetic, total, _, _ = local_energy_audit(geometry, initial['psi'], initial['g'], result)
        local.append(np.sqrt(np.sum(mesh.volumes*(magnetic-total)**2)/np.sum(mesh.volumes*(initial['magnetic_dot_induction']/MU0)**2)))
    assert errors[0]/errors[1] > 1.8 and errors[1]/errors[2] > 1.8
    assert local[0]/local[1] > 2.5 and local[1]/local[2] > 2.5
    # No axis accuracy assertion: the measured nonvanishing first-annulus
    # defect is retained in HALL_WORK_VERIFICATION and prevents integration.


def test_invalid_carriers_and_boundary_data_fail_before_mutation():
    mesh = MHDMesh(1., 1., 8, 8)
    geometry = HallCurl(mesh)
    r, _ = mesh.corners
    psi, g = .5*r*r, np.zeros((8, 8))
    saved = psi.copy()
    for density in (0., -1., np.nan, np.inf):
        with pytest.raises(ValueError):
            hall_induction(geometry, psi, g, density)
    psi[-1, 2] *= 1.01
    with pytest.raises(ValueError):
        hall_induction(geometry, psi, g, 1e20)
    psi[-1, 2] = saved[-1, 2]
    assert_allclose(psi, saved, atol=0.)
