"""Moving ideal species limits, compression and atomic-state contracts."""
import numpy as np
import pytest
from scipy.constants import atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import MHDParameters
from plasma_bench.mhd_state import NonphysicalState
from plasma_bench.toroidal_state import from_primitive
from plasma_bench.toroidal_mhd import ToroidalMHD
from plasma_bench.warm_mhd import AdiabaticSpeciesMHD
from plasma_bench.two_temperature import SpeciesThermodynamics


def fixture(order=2, scale=1e20):
    mesh = MHDMesh(.035, .12, 6, 12)
    r, z = mesh.centers
    rc, _ = mesh.corners
    mi = 39.948*atomic_mass
    rho = mi*1e20*(1+.13*np.cos(np.pi*z/.12))
    state = from_primitive(mesh, .01*rc**2/2, rho, 110., vr=200*r*(1-(r/.035)**2), vz=300*np.sin(np.pi*z/.12), magnetic_quadrature='compatible')
    thermal = SpeciesThermodynamics(mesh, mi, density_scale=scale)
    n = thermal.density(state)
    pe = 70*(n/1e20)**(5/3)
    solver = AdiabaticSpeciesMHD(mesh, state, pe, mi, MHDParameters(spatial_order=order), density_scale=scale)
    return solver


@pytest.mark.parametrize('order', [1, 2])
def test_moving_mass_species_compression_and_identical_ideal_fluid(order):
    solver = fixture(order)
    cold = ToroidalMHD(solver.mesh, solver.state, solver.parameters)
    thermal = solver.thermal
    initial_w = solver.electron_entropy/thermal.density(solver.state)
    a, du, _, _ = solver.species_rhs(solver.state, solver.electron_entropy)
    n = thermal.density(solver.state)
    dn = thermal.mapping.density(a.conserved[..., 0])/thermal.ion_mass
    pe, _ = thermal.validate(solver.state, solver.electron_entropy)
    np.testing.assert_allclose(thermal.pressure_rate(solver.state, solver.electron_entropy, a.conserved[..., 0], du),
                               solver.parameters.gamma*pe*dn/n, atol=3e-8, rtol=2e-12)
    for _ in range(8):
        dt = .5*solver.timestep_limit()
        cold.step(dt); solver.step(dt)
        np.testing.assert_array_equal(solver.state.conserved, cold.state.conserved)
        np.testing.assert_array_equal(solver.state.psi, cold.state.psi)
    np.testing.assert_allclose(solver.electron_entropy/thermal.density(solver.state), initial_w, rtol=1e-15)
    d = solver.diagnostics()
    assert abs(d['total_energy_j_change']) < 1e-16
    assert abs(d['electron_internal_energy_j']+d['ion_internal_energy_j']-d['internal_energy_j']) < 5e-17


def test_density_scale_is_a_coordinate_choice():
    a, b = fixture(scale=1e20), fixture(scale=7e18)
    for _ in range(5):
        dt = .5*a.timestep_limit()
        a.step(dt); b.step(dt)
    np.testing.assert_array_equal(a.state.conserved, b.state.conserved)
    for pa, pb in zip(a.thermal.validate(a.state, a.electron_entropy), b.thermal.validate(b.state, b.electron_entropy)):
        np.testing.assert_allclose(pa, pb, rtol=2e-15)


def test_failed_species_trial_leaves_all_accepted_state_untouched(monkeypatch):
    solver = fixture()
    state, entropy, diagnostics = solver.state.copy(), solver.electron_entropy.copy(), solver.diagnostics()
    original = solver.thermal.validate
    count = 0
    def reject_final(*args):
        nonlocal count
        count += 1
        if count == 3:
            raise NonphysicalState('forced final species failure')
        return original(*args)
    monkeypatch.setattr(solver.thermal, 'validate', reject_final)
    with pytest.raises(NonphysicalState, match='final species failure'):
        solver.step(.5*solver.timestep_limit())
    np.testing.assert_array_equal(solver.state.conserved, state.conserved)
    np.testing.assert_array_equal(solver.state.psi, state.psi)
    np.testing.assert_array_equal(solver.electron_entropy, entropy)
    assert solver.diagnostics() == diagnostics


def test_unsupported_species_partition_and_resistivity_rejected():
    solver = fixture()
    with pytest.raises(NonphysicalState, match='partition'):
        AdiabaticSpeciesMHD(solver.mesh, solver.state, 1e5*np.ones_like(solver.electron_entropy), solver.thermal.ion_mass)
    with pytest.raises(ValueError, match='heating partition'):
        AdiabaticSpeciesMHD(solver.mesh, solver.state, np.ones_like(solver.electron_entropy), solver.thermal.ion_mass, MHDParameters(resistivity=1e-8))
