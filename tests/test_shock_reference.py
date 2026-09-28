"""Independent jumps, cylindrical balance and refinement of actual MHD kernels."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.shock_reference import CylindricalShock
from plasma_bench.shock_benchmarks import standing_shock, initial_fields
from plasma_bench.coupled_hydro import gas_primitive
from plasma_bench.mhd_state import primitive


@pytest.mark.parametrize("field", [0., np.sqrt(MU0), 2*np.sqrt(MU0)])
def test_shock_has_independent_conserved_fluxes_and_entropy_increase(field):
    ref = CylindricalShock(upstream_field=field)
    upstream = (ref.upstream_density, ref.upstream_velocity, ref.upstream_pressure, field)
    downstream = ref.downstream
    assert_allclose(ref.flux_invariants(*downstream), ref.flux_invariants(*upstream), rtol=3e-14, atol=1e-14)
    rho, v, p, b = downstream
    assert p/rho**ref.gamma > ref.upstream_pressure/ref.upstream_density**ref.gamma
    assert 1 < rho/ref.upstream_density < (ref.gamma+1)/(ref.gamma-1)
    assert v*v < ref.gamma*p/rho+b*b/(MU0*rho)
    if field == 0:
        mach2 = ref.upstream_velocity**2/(ref.gamma*ref.upstream_pressure/ref.upstream_density)
        assert_allclose(rho/ref.upstream_density, (ref.gamma+1)*mach2/(2+(ref.gamma-1)*mach2), rtol=1e-14)


@pytest.mark.parametrize("downstream,r", [(True, np.linspace(.4, .49, 7)), (False, np.linspace(.51, .8, 7))])
def test_reference_satisfies_cylindrical_momentum_and_advected_invariants(downstream, r):
    ref = CylindricalShock()
    rho, v, p, b = ref.branch(r, downstream=downstream)
    assert_allclose(r*rho*v, ref.mass_flux, rtol=1e-14)
    assert_allclose(b/rho, ref.field_ratio, rtol=1e-14)
    assert_allclose(p/rho**ref.gamma, (p/rho**ref.gamma)[0], rtol=1e-14)
    enthalpy = ref.gamma/(ref.gamma-1)*p/rho+.5*v*v+b*b/(MU0*rho)
    assert_allclose(enthalpy, ref.specific_energy, rtol=2e-14)
    h = 1e-5
    def momentum(radius):
        density, velocity, pressure, field = ref.branch(radius, downstream=downstream)
        return radius*(density*velocity**2+pressure+field*field/(2*MU0))
    # Independently differentiated conservative radial momentum equation.
    derivative = (momentum(r+h)-momentum(r-h))/(2*h)
    assert_allclose(derivative, p+b*b/(2*MU0), rtol=2e-7, atol=1e-8)
    with pytest.raises(ValueError, match="choking"):
        ref.branch(.02, downstream=downstream)


@pytest.mark.parametrize("coupled", [False, True])
def test_standing_magnetic_shock_profile_error_decreases(coupled):
    errors = []
    for cells in (24, 48):
        solver, ref = standing_shock(cells, coupled=coupled, shock_diffusion=1.)
        gas, flux = initial_fields(solver.mesh, ref)
        exact = gas_primitive(solver.mesh, gas, ref.gamma)
        exact_b = solver.mesh.cell_fields(flux)[1]
        solver.advance(.02)
        if coupled:
            q = gas_primitive(solver.mesh, solver.state.gas, ref.gamma)
            bz = solver.mesh.cell_fields(solver.mapping.fluid_flux(solver.total_flux()))[1]
        else:
            q = primitive(solver.mesh, solver.state, ref.gamma)
            bz = q[..., 5]
        mask = (solver.mesh.r[:, None] > .44)&(solver.mesh.r[:, None] < .66)&(abs(solver.mesh.z[None, :]) < .1)
        errors.append([np.mean(abs((q[..., 0]-exact[..., 0])[mask])),
                       np.mean(abs((bz-exact_b)[mask]))/ref.upstream_field])
        assert bz[mask].min() > 0  # Analytic field has no reversal.
        d = solver.diagnostics()
        energy_key = "response_energy_balance_relative" if coupled else "total_energy_relative_change"
        assert abs(d[energy_key]) < 2e-12
        assert abs(d["mass_relative_change"]) < 1e-13
    assert np.all(np.asarray(errors[1]) < .8*np.asarray(errors[0]))
