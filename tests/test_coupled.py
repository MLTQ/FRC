import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.geometry import Grid
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.coupled import CoupledState, CoupledParameters, VacuumCoupledMHD
from plasma_bench.coupled_hydro import gas_from_primitive, gas_primitive, hydro_step
from plasma_bench.coupled_magnetic import MagneticStepper
from plasma_bench.coupled_benchmarks import (aligned_meshes, compact_flux, compact_force,
    cell_average, passive_advection, passive_diffusion, gaussian_diffusion, coil_pulse)
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_benchmarks import sod_exact
from plasma_bench.vacuum_benchmarks import infinite_flux, finite_box_flux


def test_mapping_has_no_exterior_fluid_and_preserves_heat():
    mesh, grid = aligned_meshes(8)
    mapping = CouplingMesh(mesh, grid)
    assert mapping.interpolate_r.shape == (mapping.node_count, mesh.nr*mesh.nz)
    assert mapping.vacuum.sum() > 0
    assert mapping.interpolate_r[mapping.vacuum].nnz == 0
    assert mapping.interpolate_z[mapping.vacuum].nnz == 0
    assert_allclose(np.asarray(mapping.heat_to_cells.sum(axis=0)).ravel(), mapping.active)
    assert mapping.heat_to_cells.data.min() > 0
    with pytest.raises(ValueError, match="two exterior"):
        CouplingMesh(mesh, Grid(1+mesh.dr, 1+mesh.dz, mesh.nr+2, mesh.nz+3))


def test_uniform_bias_and_resting_gas_are_stationary():
    mesh, grid = aligned_meshes(8)
    state = CoupledState(gas_from_primitive(mesh, 1., 1.), np.zeros(grid.shape))
    solver = VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=np.sqrt(MU0)), state)
    solver.advance(.05)
    assert_allclose(solver.state.gas, state.gas, rtol=0, atol=1e-14)
    assert_allclose(solver.state.induced_psi, 0, atol=1e-18)
    assert abs(solver.diagnostics()["response_energy_balance_relative"]) < 1e-14


def test_closed_end_gas_shock_against_independent_sod_solution():
    errors = []
    for nz in (80, 160):
        mesh = MHDMesh(1, 1, 4, nz)
        rho, pressure = np.where(mesh.z < 0, 1., .125), np.where(mesh.z < 0, 1., .1)
        state = gas_from_primitive(mesh, rho[None, :], pressure[None, :], gamma=1.4)
        initial_energy = np.sum(state[..., 3]*mesh.volumes)
        for _ in range(nz):
            state = hydro_step(mesh, state, .12/nz, 1.4)
        q = gas_primitive(mesh, state, 1.4)
        exact = sod_exact(mesh.z, .12)[0]
        select = abs(mesh.z) < .4
        errors.append(np.mean(abs(q[1, select, 0]-exact[select])))
        assert_allclose(np.sum(state[..., 3]*mesh.volumes), initial_energy, rtol=2e-14)
    assert errors[1] < .7*errors[0]
    assert errors[1] < .015


def test_independent_continuum_magnetic_force_converges():
    errors = []
    for cells in (8, 16, 32):
        mesh, grid = aligned_meshes(cells)
        mapping = CouplingMesh(mesh, grid)
        u = compact_flux(grid.r[:, None], grid.z[None, :])
        total = u+.5*np.sqrt(MU0)*grid.r[:, None]**2
        force = MagneticStepper(mapping, 0).advection(total).T@((2*np.pi/MU0)*(mapping.k@mapping.interior(u)))
        measured = [f.reshape(mesh.nr, mesh.nz)/mesh.volumes for f in np.split(force, 2)]
        exact = cell_average(mesh, compact_force)
        error = sum((a-b)**2 for a, b in zip(measured, exact))
        norm = sum(a*a for a in exact)
        errors.append(np.sqrt(np.sum(error*mesh.volumes)/np.sum(norm*mesh.volumes)))
    assert errors[1] < errors[0]/3
    assert errors[2] < errors[1]/3.4


@pytest.mark.parametrize("resistivity", [0., .01*MU0])
def test_magnetic_coil_work_lorentz_exchange_and_ohmic_heat(resistivity):
    mesh, grid = aligned_meshes(8)
    mapping = CouplingMesh(mesh, grid)
    u = compact_flux(grid.r[:, None], grid.z[None, :], amplitude=.1)
    gas = gas_from_primitive(mesh, 1., 1.)
    ext0 = .5*np.sqrt(MU0)*grid.r[:, None]**2+np.zeros((1, grid.nz))
    update = MagneticStepper(mapping, resistivity).step(gas, u, ext0, 1.02*ext0, .005)
    old_u, new_u = mapping.interior(u), mapping.interior(update.induced)
    delta_field = np.pi/MU0*(new_u@(mapping.k@new_u)-old_u@(mapping.k@old_u))
    delta_gas = np.sum((update.gas[..., 3]-gas[..., 3])*mesh.volumes)
    assert_allclose(delta_field+delta_gas, update.source_work, rtol=2e-11, atol=1e-14)
    assert update.lorentz_work > 0
    delta_internal = np.sum((gas_primitive(mesh, update.gas, 5/3)[..., 3]-1.)*mesh.volumes/(5/3-1))
    assert_allclose(delta_internal, update.ohmic_heat, atol=3e-15, rtol=1e-10)
    assert (update.ohmic_heat > 0) == (resistivity > 0)
    ku = mapping.k@new_u
    assert np.linalg.norm(ku[mapping.vacuum]) < 1e-11*np.linalg.norm(ku)


def test_passive_magnetic_advection_is_second_order():
    errors = []
    for cells in (12, 24):
        solver = passive_advection(cells)
        solver.advance(.05)
        grid = solver.mapping.field
        exact = compact_flux(grid.r[:, None], grid.z[None, :], width=.5, amplitude=1e-6, shift=.025)
        errors.append(np.linalg.norm(solver.state.induced_psi-exact)/np.linalg.norm(exact))
    assert errors[1] < errors[0]/3.5


def test_resistive_vacuum_coupling_matches_gaussian_diffusion():
    errors = []
    for cells in (12, 24):
        solver = passive_diffusion(cells)
        solver.advance(.03)
        grid = solver.mapping.field
        exact = gaussian_diffusion(grid.r[:, None], grid.z[None, :], .03)
        select = (grid.r[:, None] < .7)&(abs(grid.z[None, :]) < .7)
        errors.append(np.linalg.norm((solver.state.induced_psi-exact)[select])/np.linalg.norm(exact[select]))
    assert errors[1] < errors[0]/3.5


def test_dynamic_pulse_accounts_for_source_work_without_exterior_current():
    solver = coil_pulse(8)
    solver.advance(.1)
    d = solver.diagnostics()
    assert d["kinetic_energy_j"] > 1e-3
    assert d["ohmic_heat_j"] > .05
    assert abs(d["mass_relative_change"]) < 2e-13
    assert abs(d["response_energy_balance_relative"]) < 2e-13
    assert d["max_vacuum_induced_current_residual_a_m2"] < 1e-7
    assert d["max_face_divergence_t_m"] < 1e-12


def test_vacuum_domain_error_is_distinct_from_discretization():
    exact = infinite_flux()
    assert_allclose(infinite_flux(order=48), exact, rtol=2e-12)
    errors = []
    for cells in (8, 16):
        measured, residual = finite_box_flux(cells, exact_boundary=True)
        errors.append(np.linalg.norm(measured-exact)/np.linalg.norm(exact))
        assert residual < 1e-12
    assert errors[1] < errors[0]/3
    small, _ = finite_box_flux(16, 2)
    large, _ = finite_box_flux(16, 8)
    assert np.linalg.norm(large-exact) < np.linalg.norm(small-exact)/15


def test_failed_magnetic_iteration_leaves_accepted_state_unchanged():
    mesh, grid = aligned_meshes(8)
    u = compact_flux(grid.r[:, None], grid.z[None, :], amplitude=1.)
    state = CoupledState(gas_from_primitive(mesh, 1., 1.), u)
    params = CoupledParameters(nonlinear_tolerance=1e-15, max_nonlinear_iterations=2)
    solver = VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=np.sqrt(MU0)), state, params)
    with pytest.raises(RuntimeError, match="did not converge"):
        solver.step()
    assert_allclose(solver.state.gas, state.gas, rtol=0, atol=0)
    assert_allclose(solver.state.induced_psi, state.induced_psi, rtol=0, atol=0)
    assert solver.time == solver.steps == solver.source_work == 0


def test_source_time_cap_and_internal_filament_rejection():
    from plasma_bench.coils import Coil
    from plasma_bench.waveforms import DC
    mesh, grid = aligned_meshes(8)
    initial = CoupledState(gas_from_primitive(mesh, 1., 1.), np.zeros(grid.shape))
    parameters = CoupledParameters(max_timestep=1e-4)
    solver = VacuumCoupledMHD(mesh, VacuumField(grid, ()), initial, parameters)
    assert solver.timestep_limit() == 1e-4
    with pytest.raises(ValueError, match="outside"):
        VacuumCoupledMHD(mesh, VacuumField(grid, (Coil(.71, .31, DC(1)),)), initial)
    with pytest.raises(ValueError, match="time cap"):
        CoupledParameters(max_timestep=0.)
