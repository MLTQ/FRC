"""Independent manufactured nonautonomous ODEs exercise coupled SSP stages."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.conservative_coupled import ConservativeParameters
from plasma_bench.coupled_benchmarks import coil_pulse
from plasma_bench.flux_limiting import LimiterReport
from plasma_bench.mhd_state import NonphysicalState


@pytest.mark.parametrize("order,minimum_ratio", [(2, 3.8), (3, 7.5)])
def test_nonautonomous_stage_times_and_positive_budget_quadrature(order, minimum_ratio):
    errors = []
    for dt in (.01, .005, .0025):
        solver = from_split_initial(coil_pulse(8, vacuum_boundary="multipole"), time_order=order,
                                    energy_quadrature="staggered", boundary_poynting="variational")
        initial_energy = solver.state.conserved[..., 3].copy()
        initial_flux = solver.mapping.fluid_flux(solver.total_flux()).copy()
        empty = np.zeros((solver.mesh.nr, solver.mesh.nz), dtype=bool)
        # Manufactured E'=3E+t, psi'=-3psi+t psi(0). These are integrator
        # verification equations, not claimed MHD heating or induction laws.
        # Actual time-varying external coils and vacuum completion remain live.
        def manufactured(state, time, total, step):
            du = np.zeros_like(state.conserved)
            du[..., 3] = 3*state.conserved[..., 3]+time
            dpsi = -3*total+time*initial_flux
            power = float(np.sum(du[..., 3]*solver.mesh.volumes))
            return (du, dpsi, -power, power, 2*power, 3*power), LimiterReport(empty, empty)
        solver._forward_euler = manufactured
        solver.timestep_limit = lambda: dt
        solver.advance(.08)
        t = solver.time
        exact_energy = (initial_energy+1/9)*np.exp(3*t)-t/3-1/9
        exact_flux = ((1+1/9)*np.exp(-3*t)+t/3-1/9)*initial_flux
        energy_error = np.max(abs(solver.state.conserved[..., 3]-exact_energy))
        flux_error = np.max(abs(solver.mapping.fluid_flux(solver.total_flux())-exact_flux))/np.max(abs(initial_flux))
        errors.append([energy_error, flux_error])
        gain = np.sum((solver.state.conserved[..., 3]-initial_energy)*solver.mesh.volumes)
        assert_allclose([solver.source_work, solver.outward_energy, solver.ohmic_heat,
                         solver.numerical_interface_heat], [gain, -gain, 2*gain, 3*gain], rtol=2e-14)
        assert solver.euler_stages == order*solver.steps
        assert solver.limited_euler_stages == 0
    assert np.all(np.array(errors[:-1])/errors[1:] > minimum_ratio)


def test_third_stage_failure_commits_no_state_work_or_counters():
    solver = from_split_initial(coil_pulse(8, vacuum_boundary="multipole"), time_order=3,
                                energy_quadrature="staggered", boundary_poynting="variational")
    original = solver._forward_euler
    calls = []
    def failing(*args):
        calls.append(args[1])
        if len(calls) == 3:
            raise NonphysicalState("Manufactured third-stage rejection")
        return original(*args)
    solver._forward_euler = failing
    old, before = solver.state.copy(), solver.diagnostics()
    with pytest.raises(NonphysicalState, match="third-stage rejection"):
        solver.step()
    assert len(calls) == 3
    assert solver.diagnostics() == before
    assert_allclose(solver.state.conserved, old.conserved, rtol=0, atol=0)
    assert_allclose(solver.state.induced_psi, old.induced_psi, rtol=0, atol=0)
    assert not solver.limiter_cell_counts.any() and not solver.limiter_troubled_counts.any()


@pytest.mark.parametrize("order", [True, 1, 4])
def test_invalid_time_order_is_rejected(order):
    with pytest.raises(ValueError, match="Time order"):
        ConservativeParameters(time_order=order)
