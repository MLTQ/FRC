"""Independent exchange, stencil propagation and magnetic admissibility checks."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import from_primitive, primitive, MHDState, NonphysicalState
from plasma_bench.flux_update import FluxUpdate
from plasma_bench.flux_limiting import admissible_update, invalid_cells


def zero_flux(mesh):
    return FluxUpdate(np.zeros((mesh.nr+1, mesh.nz, 4)),
                      np.zeros((mesh.nr, mesh.nz+1, 4)),
                      np.zeros((mesh.nr+1, mesh.nz+1)), np.zeros((mesh.nr, mesh.nz)))


def uniform_state(mesh):
    r, z = mesh.corners
    return from_primitive(mesh, .5*np.sqrt(MU0)*r*r, 1., 2/3)


def test_local_replacement_rechecks_newly_invalid_neighbors_and_conserves_exchange():
    mesh = MHDMesh(1., 1., 10, 12)
    state = uniform_state(mesh)
    old = state.copy()
    dt = .1
    high, low = zero_flux(mesh), zero_flux(mesh)
    # An energy conveyor has only one invalid cell initially. Removing its
    # outflow exposes the next deficit; a single pass is insufficient.
    high.axial[3, 3:7, 3] = 2*mesh.dz/dt
    high.axial[3, 3:7, 0] = .2*mesh.dz/dt
    r, z = mesh.corners
    high.electric[:] = .001*np.sqrt(MU0)*np.sin(np.pi*r)*np.cos(np.pi*z/2)
    initial_rhs = high.rhs(mesh)
    initial_bad = invalid_cells(mesh, state.conserved+dt*initial_rhs[0], state.psi+dt*initial_rhs[1])
    assert np.count_nonzero(initial_bad) == 1 and initial_bad[3, 2]
    rhs, report = admissible_update(mesh, state.conserved, state.psi, dt, high, lambda: low)
    result = MHDState(state.conserved+dt*rhs[0], state.psi+dt*rhs[1])
    primitive(mesh, result, 5/3)
    assert report.iterations == 4
    assert np.count_nonzero(report.troubled) == 4
    assert np.count_nonzero(report.affected) < mesh.nr*mesh.nz/4
    assert_allclose(rhs[0][~report.affected], initial_rhs[0][~report.affected], rtol=0, atol=0)
    for actual, exact in zip(mesh.cell_fields(rhs[1]), mesh.cell_fields(initial_rhs[1])):
        assert_allclose(actual[~report.affected], exact[~report.affected], rtol=0, atol=0)
    assert_allclose(np.sum((result.conserved-old.conserved)*mesh.volumes[..., None], axis=(0, 1))[[0, 3]],
                    [0., -dt*rhs[2]], atol=3e-15)
    assert abs(mesh.divergence(result.psi)).max() < 3e-17
    assert_allclose(state.conserved, old.conserved, rtol=0, atol=0)
    assert_allclose(state.psi, old.psi, rtol=0, atol=0)


def test_magnetic_energy_alone_triggers_shared_corner_fallback():
    mesh = MHDMesh(1., 1., 10, 12)
    state = uniform_state(mesh)
    high, low = zero_flux(mesh), zero_flux(mesh)
    # Fluid total energy is unchanged, but a concentrated electric curl would
    # make magnetic energy exceed it. Checking U alone would miss this.
    high.electric[5, 6] = np.sqrt(MU0)
    rhs, report = admissible_update(mesh, state.conserved, state.psi, 1., high, lambda: low)
    assert 0 < report.troubled.sum() <= 4
    assert report.affected.sum() < mesh.nr*mesh.nz/4
    result = MHDState(state.conserved+rhs[0], state.psi+rhs[1])
    primitive(mesh, result, 5/3)
    assert_allclose(result.conserved, state.conserved, rtol=0, atol=0)
    assert_allclose(result.psi, state.psi, rtol=0, atol=0)


def test_valid_proposal_does_not_evaluate_baseline_or_change_fluxes():
    mesh = MHDMesh(1., 1., 8, 8)
    state = uniform_state(mesh)
    high = zero_flux(mesh)
    def forbidden():
        pytest.fail("An admissible high-order proposal must not evaluate low order")
    rhs, report = admissible_update(mesh, state.conserved, state.psi, .1, high, forbidden)
    assert report.minimum_weight == 1 and not report.affected.any()
    assert_allclose(rhs[0], 0)


def test_radial_one_dimensional_frozen_ratio_survives_local_flux_choices():
    mesh = MHDMesh(1., 1., 12, 8)
    state = uniform_state(mesh)
    high, low = zero_flux(mesh), zero_flux(mesh)
    high.radial[..., 0] = .1*np.sin(np.pi*mesh.rf[:, None])
    high.electric[:] = np.sqrt(MU0)*high.radial[:, :1, 0]
    high.radial_source[4, :] = 100.
    rhs, report = admissible_update(mesh, state.conserved, state.psi, .1, high, lambda: low)
    result = MHDState(state.conserved+.1*rhs[0], state.psi+.1*rhs[1])
    q = primitive(mesh, result, 5/3)
    assert report.troubled[4].all() and report.affected.sum() < mesh.nr*mesh.nz/2
    assert_allclose(q[..., 5]/np.sqrt(MU0), q[..., 0], atol=3e-15, rtol=0)


def test_selected_boundary_flux_controls_energy_budget():
    mesh = MHDMesh(1., 1., 8, 12)
    state = uniform_state(mesh)
    high, low = zero_flux(mesh), zero_flux(mesh)
    high.radial[-1, :, 3] = .03
    low.radial[-1, :, 3] = .01
    high.radial_source[-1, 4] = 100.
    rhs, report = admissible_update(mesh, state.conserved, state.psi, .1, high, lambda: low)
    assert report.troubled.sum() == 1 and report.troubled[-1, 4]
    expected = 2*np.pi*mesh.radius*mesh.dz*(.03*(mesh.nz-1)+.01)
    assert_allclose(rhs[2], expected, rtol=2e-15)
    assert_allclose(np.sum(rhs[0][..., 3]*mesh.volumes), -expected, rtol=2e-15)


@pytest.mark.parametrize("mode", ["local", "global"])
def test_invalid_baseline_raises_without_state_repair(mode):
    mesh = MHDMesh(1., 1., 8, 8)
    state = uniform_state(mesh)
    old = state.copy()
    bad = zero_flux(mesh)
    bad.radial_source[3, 3] = 100.
    with pytest.raises(NonphysicalState, match="First-order"):
        admissible_update(mesh, state.conserved, state.psi, 1., bad, lambda: bad, mode=mode)
    assert_allclose(state.conserved, old.conserved, rtol=0, atol=0)
    assert_allclose(state.psi, old.psi, rtol=0, atol=0)
