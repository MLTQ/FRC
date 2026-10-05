"""Physical wall traces, original variational work and actual Gauss equations."""
import numpy as np
import pytest
from examples.strong_wall_reference import initialize
from examples.current_transport_sampling import bases
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.strong_wall_coordinates import StrongWallCoordinates
from plasma_bench.strong_wall_gauss import StrongWallGauss


def fixture(kind='contact'):
    return initialize(dict(kind=kind, elements=1, points=8, field_nr=6, field_nz=8, amplitude=1e-4))


def test_independent_wall_traces_ampere_gradient_and_zero_reaction_work():
    ref, d, u, s, b = fixture(); geo = s.geometry; old = geo.original
    c = StrongWallCoordinates(s, u, b); rng = np.random.default_rng(190)
    du, db = c.unpack(1e-4*rng.normal(size=c.size)); u = u+du; b = c.field(db)
    model = RegularCurrent(s.mapping, u[..., 0], ref.ion_mass)
    original = RegularCurrent(s.original_system.mapping, u[..., 0], ref.ion_mass)
    j = model.current(db); j0 = original.current(db)
    assert np.linalg.norm(j-j0)/np.linalg.norm(j0) < 3e-14
    # Independent root-product physical boundary interpolation, not w@P.
    trace = bases(.5*geo.mesh.rf**2, np.array([.5*d.radius**2]), 5)[3]
    jr = j[:geo.mesh.nr*geo.mesh.nz].reshape(geo.mesh.nr, geo.mesh.nz)
    assert np.linalg.norm(trace@jr)/np.linalg.norm(jr) < 2e-12
    assert np.linalg.norm(trace@b.g)/np.linalg.norm(b.g) < 2e-12
    assert np.all(u[-1, -1, :, 1] == 0)
    potential = rng.normal(size=(geo.mesh.nr+1, geo.mesh.nz))
    grad = geo.gradient(potential)
    assert np.linalg.norm(geo.faraday@grad)/(np.linalg.norm(geo.faraday.toarray())*np.linalg.norm(grad)) < 2e-16
    rate = s.derivative(u, b)
    raw_defect = old.faraday@rate.electric-rate.magnetic
    _, test = c.unpack(rng.normal(size=c.size))
    work = test@(geo.hessian@raw_defect)
    scale = np.sqrt(test@(geo.hessian@test))*np.sqrt(raw_defect@(geo.hessian@raw_defect))
    assert abs(work)/scale < 1e-9
    # Independently differentiate physical grid+quadrature energy, including n_t.
    fluid_power = np.sum(d.volume*np.sum(d.fields(u)['energy_gradient']*rate.fluid, axis=-1))
    current_power = j@(model.inertia@rate.current)+np.sum(model.density_energy_gradient(j)*rate.fluid[..., 0])
    magnetic_power = db@(geo.hessian@rate.magnetic)
    scale = abs(fluid_power)+abs(current_power)+abs(magnetic_power)
    assert abs(fluid_power+current_power+magnetic_power)/scale < 1e-10
    assert rate.residual <= s.solver_options['tolerance']
    assert np.all(rate.fluid[-1, -1, :, 1] == 0)


def test_full_wall_gauss_equations_and_atomic_rejection():
    ref, _, u, s, b = fixture('wave'); solver = StrongWallGauss(s, u, b)
    c = solver.coordinates; before = u.copy(); initial = b.packed_increment(s.geometry)
    dt = 2*np.pi/abs(ref.omega)/60
    v, field, diagnostic = solver.step(u, b, dt)
    rates = []; actual = []
    for y in diagnostic['stage_changes']:
        du, db = c.unpack(y); state = u+du; magnetic = initial+db
        actual.append(c.pack(state-u, magnetic-initial))
        stage_field = c.field(magnetic)
        assert np.all(state[-1, -1, :, 1] == 0)
        assert s.geometry.wall_residual(stage_field.g) < 256*np.finfo(float).eps
        rate = s.derivative(state, stage_field); rates.append(c.pack(rate.fluid, rate.magnetic))
    residual = np.asarray(actual)-dt*solver.a@rates
    endpoint = c.pack(v-u, field.packed_increment(s.geometry)-initial)-.5*dt*np.sum(rates, axis=0)
    assert max(np.max(abs(residual)), np.max(abs(endpoint))) <= diagnostic['residual_limit']
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(s.geometry), initial)
    ref, _, u, s, b = fixture(); before = u.copy(); initial = b.packed_increment(s.geometry)
    exhausted = StrongWallGauss(s, u, b, max_iterations=1)
    with pytest.raises(RuntimeError, match='did not converge'):
        exhausted.step(u, b, ref.period/120)
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(s.geometry), initial)
    bad = u.copy(); bad[-1, -1, 0, 1] = 1.
    with pytest.raises(ValueError, match='wall-normal'): s.derivative(bad, b)
