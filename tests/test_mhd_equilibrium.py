import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_equilibrium import PeriodicEquilibrium
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_state import NonphysicalState


@pytest.mark.parametrize("aspect", [.5, 1., 2., 4.])
def test_independent_equilibrium_series_and_continuum_force_balance(aspect):
    eq = PeriodicEquilibrium(radius=.7, half_length=.7*aspect)
    assert abs(eq.radial_flux(1)) < 1e-12
    assert max(abs(eq.ode_residual(np.linspace(0, 1, 301)))) < 1e-10
    r, z = np.meshgrid(np.linspace(.12, .61, 11), np.linspace(-.23, .23, 13), indexing="ij")
    h = 2e-6
    _, br, bz, _ = eq.fields(r, z)
    rp, rm, zp, zm = eq.fields(r+h, z), eq.fields(r-h, z), eq.fields(r, z+h), eq.fields(r, z-h)
    current = ((zp[1]-zm[1])-(rp[2]-rm[2]))/(2*h*MU0)
    grad_p_r, grad_p_z = (rp[3]-rm[3])/(2*h), (zp[3]-zm[3])/(2*h)
    force = np.stack((current*bz, -current*br))
    gradient = np.stack((grad_p_r, grad_p_z))
    assert np.linalg.norm(force-gradient)/np.linalg.norm(gradient) < 2e-8
    assert eq.fields(0, 0)[2] > 0
    assert eq.fields(eq.radius, 0)[2] < 0


def test_curved_equilibrium_drift_and_force_decrease_with_refinement():
    # The default background-pressure case is resolved on inexpensive meshes.
    # A lower-pressure stress fixture's coarse-grid failures are separately
    # recorded and must never be converted into floors or hidden repairs.
    eq = PeriodicEquilibrium()
    forces, velocities, flux_errors = [], [], []
    for n in (24, 48):
        mesh = MHDMesh(1, 1, n, 2*n)
        initial = eq.state(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        rhs, _ = solver.rhs(initial)
        forces.append(np.sqrt(np.sum(np.sum(rhs.conserved[..., 1:3]**2, axis=-1)*mesh.volumes)/mesh.volumes.sum()))
        solver.advance(.05)
        velocities.append(solver.diagnostics()["max_velocity_m_s"])
        flux_errors.append(np.linalg.norm(solver.state.psi-initial.psi)/np.linalg.norm(initial.psi))
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
        assert abs(solver.diagnostics()["mass_relative_change"]) < 1e-12
        field = max(abs(v).max() for v in mesh.face_fields(solver.state.psi))
        assert np.max(abs(mesh.divergence(solver.state.psi)))*mesh.dr/field < 1e-13
    assert forces[1] < .4*forces[0]
    assert velocities[1] < .5*velocities[0]
    assert flux_errors[1] < .4*flux_errors[0]


def test_unresolved_equilibrium_reconstruction_fails_without_state_repair():
    mesh = MHDMesh(1, 1, 12, 24)
    initial = PeriodicEquilibrium(background_pressure=.1).state(mesh)
    solver = AxisymmetricMHD(mesh, initial, MHDParameters(positivity_reconstruction=False))
    with pytest.raises(NonphysicalState, match="Reconstructed"):
        solver.step()
    assert_allclose(solver.state.conserved, initial.conserved, atol=0, rtol=0)
    assert_allclose(solver.state.psi, initial.psi, atol=0, rtol=0)
    assert solver.time == 0
