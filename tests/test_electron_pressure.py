"""Independent current/pressure/thermal conservation and evolution checks."""
import numpy as np
import pytest
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.integrate import solve_ivp
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_current import DualCurrent
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_midpoint import step, ElectronConvergenceError
from plasma_bench.hall_compatible import CompatibleHallCurl
from plasma_bench.thermal_means import logarithmic_mean, entropy_mean
from examples.electron_reference import ElectronReference, drift_pressure


def fixture(nr=12):
    ref = ElectronReference()
    mesh = MHDMesh(ref.radius, ref.length/2, nr, 2*nr)
    geo = DualCurrent(mesh)
    data = ref.fields(geo.r, geo.z)
    model = ElectronPressure(mesh, data['density'])
    return ref, mesh, model, ref.cell_average(mesh, 'g'), model.from_pressure(data['pressure'])


def test_dual_current_exact_local_charge_and_gradient_work():
    mesh = MHDMesh(.035, .12, 9, 14)
    geo = DualCurrent(mesh)
    rng = np.random.default_rng(31)
    g = rng.normal(size=(mesh.nr, mesh.nz))*.04
    wall = rng.normal(size=mesh.nz)*.001
    ir, iz = geo.current_flux(g, wall)
    residual = np.diff(ir, axis=0)+iz-np.roll(iz, 1, axis=1)
    assert abs(residual).max() < 2e-15*max(abs(ir).max(), abs(iz).max())
    phi = rng.normal(size=(mesh.nr+1, mesh.nz))
    er, ez = geo.gradient(phi)
    assert abs(geo.faraday(er, ez)).max() < 2e-15*max(abs(er).max()/mesh.dz, abs(ez).max()/geo.dx.min())
    er, ez = rng.normal(size=er.shape), rng.normal(size=ez.shape)
    work = geo.work(g, wall, er, ez)
    assert abs(work['magnetic_rate_w']+work['current_work_w']-work['boundary_magnetic_input_w']) < 1e-13*max(abs(v) for v in work.values())
    assert np.sum(geo.volumes) == pytest.approx(np.pi*mesh.radius**2*2*mesh.half_length, rel=3e-15)


def test_dual_current_matches_hall_adjoint_and_uniform_physical_current():
    _, mesh, model, g, _ = fixture(8)
    geo = model.geometry
    hall = CompatibleHallCurl(mesh)
    psi = np.zeros((mesh.nr+1, mesh.nz+1))
    dual = hall.split_electric(-hall.faraday.T @ hall.magnetic_gradient(psi, g))
    ir, iz = geo.current_flux(g, 0.)
    np.testing.assert_allclose(dual[0].reshape(mesh.nr, mesh.nz), ir[1:-1]*geo.dx[:, None], rtol=4e-15, atol=3e-16)
    np.testing.assert_allclose(dual[2].reshape(mesh.nr, mesh.nz), iz[:-1]*mesh.dz, rtol=4e-15, atol=3e-15)
    c = .025
    ir, iz = geo.current_flux(np.full_like(g, c), c*mesh.radius)
    assert np.max(abs(ir)) == 0
    np.testing.assert_allclose(iz/(geo.volumes/mesh.dz), 2*c/MU0, rtol=7e-15)
    _, omitted = geo.current_flux(np.full_like(g, c), 0.)
    np.testing.assert_allclose(omitted[-1]-iz[-1], -2*np.pi*mesh.radius**2*c/MU0, rtol=3e-15)
    assert np.all(omitted[-1] < 0) and np.all(iz[-1] > 0)


def test_independent_battery_and_thermal_rate_refinement():
    errors, axis = [], []
    for nr in (8, 16, 32):
        ref, mesh, model, g, u = fixture(nr)
        dg, du, _, _ = model.rhs(g, u)
        exact = ref.cell_average(mesh, 'g_rate')
        exact_u = ref.fields(model.geometry.r, model.geometry.z)['entropy_rate']
        errors.append((np.sqrt(np.sum(mesh.volumes*(dg-exact)**2)/np.sum(mesh.volumes*exact**2)),
                       np.sqrt(np.sum(model.geometry.volumes*(du-exact_u)**2)/np.sum(model.geometry.volumes*exact_u**2))))
        axis.append(np.linalg.norm(dg[0]-exact[0])/np.linalg.norm(exact[0]))
    assert np.all(np.array(errors[:-1])/np.array(errors[1:]) > 3.)
    assert axis[0]/axis[1] > 3. and axis[1]/axis[2] > 3.


def test_paired_pressure_energy_and_entropy_proxy_invariants():
    _, mesh, model, g, u = fixture(12)
    dg, du, er, ez = model.rhs(g, u)
    geo = model.geometry
    np.testing.assert_allclose(dg, geo.faraday(er, ez), rtol=2e-10, atol=2e-7)
    thermal = np.sum(geo.volumes*model.energy_weight*du)
    magnetic = np.sum(geo.magnetic_mass*g*dg)
    scale = np.sum(abs(geo.volumes*model.energy_weight*du))+np.sum(abs(geo.magnetic_mass*g*dg))
    assert abs(thermal+magnetic) < 4e-15*scale
    assert abs(np.sum(geo.volumes*du)) < 1e-14*np.sum(geo.volumes*abs(du))
    assert abs(np.sum(geo.volumes*model.density/u*du)) < 3e-13*np.sum(geo.volumes*abs(model.density/u*du))
    before = model.diagnostics(g, u)
    out = step(model, g, u, 1e-8, tolerance=2e-13)
    after = model.diagnostics(out.g, out.entropy_density)
    assert abs(after['total_energy_j']/before['total_energy_j']-1) < 4e-15
    particles = np.sum(geo.volumes*model.density)
    assert abs(after['entropy_log_particle_units']-before['entropy_log_particle_units'])/particles < 3e-15
    assert np.linalg.norm(out.g-g) > 0 and np.linalg.norm(out.entropy_density-u) > 0


def test_adiabatic_barotrope_and_density_scale_invariance():
    _, mesh, model, g, u = fixture(10)
    pressure = model.density*model.a*5*E
    dg, _, er, ez = model.rhs(g, model.from_pressure(pressure))
    assert abs(dg).max() < 2e-14*max(abs(er).max()/mesh.dz, abs(ez).max()/model.geometry.dx.min())
    dg, _, er, ez = model.rhs(g, model.from_pressure(model.density*5*E))
    assert abs(dg).max() < 2e-12*max(abs(er).max()/mesh.dz, abs(ez).max()/model.geometry.dx.min())
    second = ElectronPressure(mesh, model.density, density_scale=7e18)
    pressure = model.pressure(u)
    dg1, du1, er1, ez1 = model.rhs(g, u)
    dg2, du2, er2, ez2 = second.rhs(g, second.from_pressure(pressure))
    np.testing.assert_allclose(dg1, dg2, rtol=2e-11, atol=1e-8)
    np.testing.assert_allclose(model.a*du1, second.a*du2, rtol=2e-12, atol=1e-8)


def test_uniform_current_electron_heat_translation():
    errors = []
    for nz in (16, 32, 64):
        mesh = MHDMesh(.035, .12, 4, nz)
        geo = DualCurrent(mesh)
        n, c = 1e20, .01
        model = ElectronPressure(mesh, np.full_like(geo.r, n), wall_btheta=c*mesh.radius)
        u = model.from_pressure(drift_pressure(geo, 0.))
        g = np.full((mesh.nr, mesh.nz), c)
        end = .24/(2*c/(MU0*n*E))
        for _ in range(4*nz):
            out = step(model, g, u, end/(4*nz), tolerance=2e-13)
            g, u = out.g, out.entropy_density
        exact = drift_pressure(geo, end)
        errors.append(np.linalg.norm(model.pressure(u)-exact)/np.linalg.norm(exact-80.))
        np.testing.assert_allclose(g, c, rtol=1e-11, atol=1e-13)
    assert errors[0]/errors[1] > 3.8 and errors[1]/errors[2] > 3.8


def test_midpoint_time_order_against_independent_adaptive_ode():
    _, mesh, model, g, u = fixture(4)
    count = g.size
    initial = np.r_[g.ravel(), u.ravel()]
    def rhs(_, y):
        dg, du, _, _ = model.rhs(y[:count].reshape(g.shape), y[count:].reshape(u.shape))
        return np.r_[dg.ravel(), du.ravel()]
    end = 5e-7
    ref = solve_ivp(rhs, (0, end), initial, method='DOP853', rtol=2e-12, atol=2e-14).y[:, -1]
    errors = []
    for steps in (4, 8, 16):
        a, b = g.copy(), u.copy()
        for _ in range(steps):
            out = step(model, a, b, end/steps, tolerance=2e-13)
            a, b = out.g, out.entropy_density
        errors.append(np.linalg.norm(a.ravel()-ref[:count]))
    assert errors[0]/errors[1] > 3.8 and errors[1]/errors[2] > 3.8


def test_invalid_closure_and_failed_midpoint_leave_inputs_unchanged():
    _, mesh, model, g, u = fixture(6)
    oldg, oldu = g.copy(), u.copy()
    with pytest.raises(ElectronConvergenceError):
        step(model, g, u, 1e-6, tolerance=1e-15, max_iterations=1)
    np.testing.assert_array_equal(g, oldg); np.testing.assert_array_equal(u, oldu)
    bad = model.density.copy(); bad[-1, 0] *= 1.1
    with pytest.raises(ValueError, match='wall density'):
        ElectronPressure(mesh, bad)
    bad[0, 0] = 0
    with pytest.raises(ValueError, match='Positive'):
        ElectronPressure(mesh, bad)


def test_thermodynamic_means_against_high_precision_reference():
    from decimal import Decimal, localcontext
    pairs = [(1., 1.), (1., 1.+1e-12), (2e-19, 2.0002e-19), (.01, 500.), (1e-100, 1e100)]
    for a, b in pairs:
        with localcontext() as context:
            context.prec = 70
            x, y = Decimal(a), Decimal(b)
            expected = x if x == y else (y-x)/(y.ln()-x.ln())
            harmonic = x*y/expected
        assert logarithmic_mean(a, b) == pytest.approx(float(expected), rel=2e-13, abs=0)
        assert entropy_mean(a, b) == pytest.approx(float(harmonic), rel=2e-13, abs=0)
    with pytest.raises(ValueError):
        logarithmic_mean(0, 1)
