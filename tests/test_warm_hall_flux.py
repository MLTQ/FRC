"""Independent magnetic boundary, charge and warm energy transport checks."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from scipy.special import j1, jn_zeros
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_current import DualCurrent
from plasma_bench.dual_mass import DualMass
from plasma_bench.hall_trace import TraceHallCurl
from plasma_bench.hall_midpoint import HallMidpoint
from plasma_bench.electron_pressure import ElectronPressure
from plasma_bench.electron_energy import energy_rate
from plasma_bench.dual_enthalpy import scalar_current_flux
from examples.electron_reference import ElectronReference


def test_full_gradient_curl_and_dual_charge_pairing_including_wall():
    mesh = MHDMesh(.035, .12, 7, 10)
    geo, dual = TraceHallCurl(mesh), DualCurrent(mesh)
    rng = np.random.default_rng(632)
    phi = rng.normal(size=(8, 10))
    er, ez = dual.gradient(phi)
    grad = np.r_[er.ravel(), np.zeros(geo.maps[1].shape[1]), ez.ravel()]
    assert np.linalg.norm(geo.faraday @ grad) < 8e-16*np.linalg.norm(geo.faraday.data)*np.linalg.norm(grad)
    g = rng.normal(size=(7, 10))
    psi = np.zeros((8, 11))
    jr, _, jz = geo.split_electric(-geo.faraday.T @ geo.magnetic_gradient(psi, g))
    ir, iz = dual.current_flux(g, 0.)
    np.testing.assert_allclose(jr.reshape(7, 10), ir[1:-1]*dual.dx[:, None], rtol=2e-14, atol=1e-12)
    np.testing.assert_allclose(jz.reshape(8, 10), iz*mesh.dz, rtol=2e-14, atol=1e-11)
    radial, axial = scalar_current_flux(dual, ir, iz, np.ones((8, 10)))
    error = np.diff(radial, axis=0)+np.diff(axial, axis=1)
    assert abs(error).max() < 8e-16*(abs(radial).max()+abs(axial).max())
    assert np.all(radial[[0, -1]] == 0)


def test_zero_magnetic_trace_allows_variable_wall_density_without_boundary_work():
    mesh = MHDMesh(.035, .12, 8, 16)
    geo = DualCurrent(mesh)
    n = 1e20*np.exp(.3*np.cos(np.pi*geo.z/.12))
    pressure = 80*(1+.2*np.sin(np.pi*geo.z/.12))
    g = .05*(1-(mesh.r[:, None]/mesh.radius)**2)**2*np.cos(np.pi*mesh.z[None]/.12)
    model = ElectronPressure(mesh, n)
    u = model.from_pressure(pressure)
    dg, du, er, ez = model.rhs(g, u)
    assert np.max(abs(ez[-1])) > 1
    work = geo.work(g, 0., er, ez)
    assert work['boundary_magnetic_input_w'] == 0
    energy = np.sum(geo.magnetic_mass*g*dg)+np.sum(geo.volumes*model.energy_weight*du)
    scale = np.sum(abs(geo.magnetic_mass*g*dg))+np.sum(abs(geo.volumes*model.energy_weight*du))
    assert abs(energy) < 2e-14*scale


def test_zero_trace_diffusion_matches_independent_dirichlet_bessel_mode():
    errors = []
    for nr in (16, 32, 64):
        mesh = MHDMesh(.035, .12, nr, 4)
        geo = TraceHallCurl(mesh)
        root = jn_zeros(1, 1)[0]/mesh.radius
        nodes, weights = np.polynomial.legendre.leggauss(16)
        r = mesh.r[:, None]+.5*mesh.dr*nodes
        g1 = .001*np.sum(weights*j1(root*r), axis=1)*mesh.dr/(2*mesh.radial_measure)
        g = np.broadcast_to(g1[:, None], (nr, 4)).copy()
        eta = 2e-5; end = .02*MU0/(eta*root**2)
        solver = HallMidpoint(geo, tolerance=1e-12)
        out = solver.step(np.zeros((nr+1, 5)), g, None, end, resistivity=eta)
        exact = g*np.exp(-.02)
        errors.append(np.linalg.norm(out.g-exact)/np.linalg.norm(exact))
        h0 = np.sum(mesh.volumes*geo.energy_density(np.zeros((nr+1, 5)), g))
        h1 = np.sum(mesh.volumes*geo.energy_density(out.psi, out.g))
        assert abs(h1-h0+out.physical_heat_j) < 2e-13*h0
        assert abs(np.sum(mesh.volumes*out.total_energy_change)) < 2e-15*h0
    assert errors[1] < .4*errors[0] and errors[2] < .5*errors[1]


def test_pressure_energy_flux_converges_to_independent_local_species_work():
    ref = ElectronReference()
    errors, defects = [], []
    for nr in (8, 16, 32):
        mesh = MHDMesh(ref.radius, ref.length/2, nr, 2*nr)
        mapping = DualMass(mesh); geo = mapping.geometry
        fields = ref.fields(geo.r, geo.z)
        model = ElectronPressure(mesh, fields['density'])
        u = model.from_pressure(fields['pressure']); g = ref.cell_average(mesh, 'g')
        actual = energy_rate(model, g, u)
        dg, du, _, _ = model.rhs(g, u)
        local = geo.magnetic_mass*g*dg/mesh.volumes+mapping.primal_average(model.energy_weight*du)
        x, w = np.polynomial.legendre.leggauss(10)
        expected = np.zeros_like(g)
        lo, hi = .5*mesh.rf[:-1]**2, .5*mesh.rf[1:]**2
        for a, wa in zip(x, w):
            r = np.sqrt(lo+hi+a*(hi-lo))[:, None]
            for b, wb in zip(x, w):
                f = ref.fields(r, mesh.z[None]+.5*mesh.dz*b)
                expected += .25*wa*wb*(r*r*f['g']*f['g_rate']/MU0+(f['density']/ref.density)**(ref.gamma-1)/(ref.gamma-1)*f['entropy_rate'])
        norm = lambda q: np.sqrt(np.sum(mesh.volumes*q*q)/np.sum(mesh.volumes*expected*expected))
        errors.append(norm(actual-expected)); defects.append(norm(actual-local))
        assert abs(np.sum(mesh.volumes*actual)) < 2e-14*np.sum(mesh.volumes*abs(actual))
    assert errors[1] < .3*errors[0] and errors[2] < .3*errors[1]
    assert defects[1] < .3*defects[0] and defects[2] < .3*defects[1]
