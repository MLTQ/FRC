"""Independent physical and boundary checks for unsplit vacuum-coupled MHD."""
from dataclasses import replace
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.coupled_benchmarks import (aligned_meshes, passive_advection, passive_diffusion,
                                            compact_flux, gaussian_diffusion, coil_pulse)
from plasma_bench.coupled_hydro import gas_from_primitive
from plasma_bench.coupled import CoupledState, VacuumCoupledMHD, CoupledParameters
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.shock_benchmarks import standing_shock, initial_fields
from plasma_bench.vacuum_trace import tangential_field, interface_resistivity
from plasma_bench.mhd_state import NonphysicalState
from plasma_bench.vacuum_flux import flux_rhs


@pytest.mark.parametrize("boundary", ["fixed", "multipole"])
def test_uniform_gas_and_bias_are_stationary(boundary):
    mesh, grid = aligned_meshes(8)
    split = VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=np.sqrt(MU0)),
                            CoupledState(gas_from_primitive(mesh, 1., 1.), np.zeros(grid.shape)),
                            CoupledParameters(vacuum_boundary=boundary))
    solver = from_split_initial(split)
    initial = solver.state.copy()
    solver.advance(.03)
    assert_allclose(solver.state.conserved, initial.conserved, atol=2e-14, rtol=0)
    assert_allclose(solver.state.induced_psi, 0, atol=1e-18)
    d = solver.diagnostics()
    assert d["limited_euler_stages"] == 0
    assert abs(d["response_energy_balance_relative"]) < 1e-14


def test_vacuum_side_trace_matches_independent_harmonic_polynomial():
    errors = []
    a, b, c = .7, .3, .05
    for n in (8, 16, 32):
        mesh, grid = aligned_meshes(n)
        mapping = CouplingMesh(mesh, grid)
        r, z = grid.r[:, None], grid.z[None, :]
        psi = a*r*r+b*r*r*z+c*r*r*(r*r-4*z*z)
        bz, low, high = tangential_field(mapping, psi)
        exact_bz = 2*a+2*b*mesh.z+c*(4*mesh.radius**2-8*(mesh.z**2+mesh.dz**2/12))
        errors.append(np.max(abs(bz-exact_bz)))
        assert_allclose(low, mesh.radial_centroid*(-b-8*c*mesh.half_length), atol=2e-13)
        assert_allclose(high, mesh.radial_centroid*(-b+8*c*mesh.half_length), atol=2e-13)
    assert errors[1] < errors[0]/3.5
    assert errors[2] < errors[1]/3.5


def test_interface_diffusivity_vanishes_with_mesh_and_is_absent_from_vacuum():
    maxima = []
    for n in (8, 16):
        solver = from_split_initial(coil_pulse(n, vacuum_boundary="multipole"))
        eta = interface_resistivity(solver.mapping, solver._primitive(solver.state, 0), 5/3)
        assert np.all(eta[solver.mapping.vacuum] == 0)
        assert np.any(eta[solver.mapping.active] == 0)
        assert eta.min() >= 0
        maxima.append(eta.max())
    assert_allclose(maxima[0]/maxima[1], 2., rtol=1e-14)


def test_interface_electric_transports_energy_on_adjacent_interior_faces():
    # Independent instantaneous Poynting identity: J=0 in uniform B, so a
    # prescribed electric curl changes magnetic, not thermal, energy at t=0.
    mesh, _ = aligned_meshes(8)
    field, alpha = np.sqrt(MU0), .01*np.sqrt(MU0)
    psi = .5*field*mesh.rf[:, None]**2+np.zeros((1, mesh.nz+1))
    q = np.zeros((mesh.nr, mesh.nz, 6))
    q[..., 0] = q[..., 3] = 1.; q[..., 5] = field
    electric = np.zeros_like(psi); electric[:, 0] = alpha*mesh.rf
    boundary = (np.full(mesh.nz, field), np.zeros(mesh.nr), np.zeros(mesh.nr))
    du, dpsi, _ = flux_rhs(mesh, q, psi, 5/3, np.zeros_like(psi), boundary, electric)
    baseline_u, baseline_psi, _ = flux_rhs(mesh, q, psi, 5/3, np.zeros_like(psi), boundary, np.zeros_like(psi))
    # Isolate the manufactured nonideal term from roundoff in reconstructing
    # an otherwise uniform B by differences of quadratic nodal flux.
    du = du-baseline_u
    dbz = mesh.cell_fields(dpsi-baseline_psi)[1]
    assert_allclose(du[..., 3]-field*dbz/MU0, 0, atol=2e-15)
    assert_allclose(du[:, 0, 3], -alpha*field/MU0, atol=2e-15)
    assert_allclose(du[:, 1:, 3], 0, atol=2e-15)


@pytest.mark.parametrize("kind", ["advection", "diffusion"])
def test_independent_smooth_solutions_converge_without_flux_limiting(kind):
    errors = []
    for n in (12, 24):
        solver = from_split_initial(passive_advection(n) if kind == "advection" else passive_diffusion(n))
        solver.advance(.05 if kind == "advection" else .03)
        grid = solver.mapping.field
        r, z = grid.r[:, None], grid.z[None, :]
        exact = compact_flux(r, z, width=.5, amplitude=1e-6, shift=.025) if kind == "advection" else gaussian_diffusion(r, z, .03)
        select = (r < .7)&(abs(z) < .7)
        errors.append(np.linalg.norm((solver.state.induced_psi-exact)[select])/np.linalg.norm(exact[select]))
        assert solver.limited_euler_stages == 0
        assert abs(solver.diagnostics()["mass_relative_change"]) < 2e-14
    assert errors[1] < errors[0]/3.2


def test_original_shock_preserves_interior_frozen_flux_and_refines_response_budget():
    errors, budgets, defects = [], [], []
    for n in (24, 48):
        split, reference = standing_shock(n)
        solver = from_split_initial(split)
        _, flux = initial_fields(solver.mesh, reference)
        exact_b = solver.mesh.cell_fields(flux)[1]
        solver.advance(.02)
        q = solver._primitive(solver.state, solver.time)
        mask = (solver.mesh.r[:, None] > .44)&(solver.mesh.r[:, None] < .66)&(abs(solver.mesh.z[None, :]) < .1)
        errors.append(np.mean(abs(q[..., 5]-exact_b)[mask])/reference.upstream_field)
        defects.append(np.mean(abs(q[..., 5]/reference.field_ratio-q[..., 0])[mask]))
        assert q[..., 5][mask].min() > 0
        d = solver.diagnostics()
        budgets.append(abs(d["response_energy_balance_relative"]))
        assert abs(d["fluid_energy_balance_residual_j"]) < 1e-11
        assert d["numerical_interface_heat_j"] > 0 and d["ohmic_heat_j"] == 0
    assert errors[1] < .8*errors[0]
    assert defects[1] < defects[0]/10 and defects[1] < 1e-6
    assert budgets[1] < .6*budgets[0]


def test_invalid_boundary_update_leaves_accepted_state_and_budgets_intact():
    solver = from_split_initial(standing_shock(24)[0], interface_diffusion=False, positivity_flux_limiting=False)
    failed = False
    while solver.time < .02:
        old, before = solver.state.copy(), solver.diagnostics()
        try:
            solver.step(min(solver.timestep_limit(), .02-solver.time))
        except NonphysicalState:
            assert_allclose(solver.state.conserved, old.conserved, rtol=0, atol=0)
            assert_allclose(solver.state.induced_psi, old.induced_psi, rtol=0, atol=0)
            assert solver.diagnostics() == before
            failed = True; break
    assert failed


def test_coil_work_error_decreases_and_local_energy_flux_is_conservative():
    errors = []
    for n in (16, 32):
        solver = from_split_initial(coil_pulse(n, vacuum_boundary="multipole"))
        solver.advance(.1)
        d = solver.diagnostics()
        errors.append(abs(d["response_residual_over_source_work"]))
        assert d["source_work_j"] > 0 and d["ohmic_heat_j"] > 0
        assert abs(d["mass_relative_change"]) < 2e-14
        assert abs(d["fluid_energy_balance_residual_j"]) < 2e-13
        assert d["max_face_divergence_t_m"] < 1e-13
    assert errors[1] < .5*errors[0]
