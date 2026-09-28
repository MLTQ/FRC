import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0
from scipy.special import j1, jn_zeros
from plasma_bench.geometry import Grid
from plasma_bench.numerics import stiffness
from plasma_bench.plasma import StationaryDiffusion, ResistiveMHD


def eigenmode(grid):
    r, z = grid.mesh
    k = jn_zeros(1, 1)[0]/grid.radius
    p = r*j1(k*r)*np.cos(np.pi*z/(2*grid.half_length))
    p[[0, -1], :] = 0
    p[:, [0, -1]] = 0
    return p, k*k+(np.pi/(2*grid.half_length))**2


def test_stiffness_symmetry_positive_energy():
    grid = Grid(nr=11, nz=21)
    k = stiffness(grid)
    assert np.max(np.abs((k-k.T).data), initial=0) < 1e-12
    u = np.random.default_rng(4).normal(size=k.shape[0])
    assert u@(k@u) > 0


def test_diffusion_eigenmode_space_and_time_convergence():
    sigma = 2e4
    errors = []
    for n in (17, 33, 65):
        g = Grid(0.05, 0.1, n, 2*n-1)
        p, eigenvalue = eigenmode(g)
        lam = eigenvalue/(mu_0*sigma)
        dt = .2/lam/100
        solver = StationaryDiffusion(g, sigma, np.zeros(g.shape), p)
        for _ in range(100):
            last = solver.step(np.zeros(g.shape), dt)
        # Compare to exact time-discrete eigenmode to isolate spatial error.
        expected = p*(1+dt*lam)**-100
        errors.append(np.linalg.norm(last.induced_psi-expected)/np.linalg.norm(expected))
    assert errors[0]/errors[1] > 3.3
    assert errors[1]/errors[2] > 3.3
    g = Grid(0.05, 0.1, 65, 129)
    p, eigenvalue = eigenmode(g)
    lam = eigenvalue/(mu_0*sigma)
    time_errors = []
    for steps in (5, 10, 20):
        solver = StationaryDiffusion(g, sigma, np.zeros(g.shape), p)
        for _ in range(steps):
            last = solver.step(np.zeros(g.shape), .5/lam/steps)
        time_errors.append(np.linalg.norm(last.induced_psi-p*np.exp(-.5))/np.linalg.norm(p*np.exp(-.5)))
    assert time_errors[0]/time_errors[1] > 1.8
    assert time_errors[1]/time_errors[2] > 1.8


def test_coupled_vacuum_ohm_ampere_and_energy_identity():
    g = Grid(0.1, .2, 31, 61)
    r, z = g.mesh
    sigma = np.where((r < .035) & (np.abs(z) < .08), 3e4, 0)
    solver = StationaryDiffusion(g, sigma, .005*r*r/2)
    u_energy = 0
    for n in range(1, 21):
        step = solver.step((.005+.01*np.sin(n/8))*r*r/2, 1e-6)
        ku = solver.k@solver.u
        ampere = ku/(mu_0*g.dr*g.dz)
        assert_allclose(ampere, step.j_theta[1:-1, 1:-1].ravel(), rtol=1e-10, atol=1e-7)
        assert np.all(step.j_theta[sigma == 0] == 0)
        assert step.ohmic_heat >= 0 and step.numerical_loss >= 0
        scale = max(abs(step.source_work), step.induced_energy, u_energy, 1e-20)
        assert abs(step.balance_residual)/scale < 2e-12
        u_energy = step.induced_energy


def test_passive_decay_and_vacuum_limit():
    g = Grid(0.05, .1, 21, 41)
    r, _ = g.mesh
    solver = StationaryDiffusion(g, 0, np.zeros(g.shape))
    assert np.all(solver.step(r*r, 1e-6).induced_psi == 0)
    p, _ = eigenmode(g)
    solver = StationaryDiffusion(g, 1e4, np.zeros(g.shape), p)
    initial = solver.induced_energy
    last = solver.step(np.zeros(g.shape), 1e-6)
    assert last.induced_energy < initial
    assert_allclose(initial-last.induced_energy, last.ohmic_heat+last.numerical_loss, rtol=1e-12)


def test_dc_produces_no_current_and_mhd_is_explicit():
    g = Grid(nr=11, nz=21)
    r, _ = g.mesh
    solver = StationaryDiffusion(g, 1e4, .005*r*r/2)
    result = solver.step(.005*r*r/2, 1e-6)
    assert np.all(result.j_theta == 0)
    with pytest.raises(NotImplementedError, match="Stage 2"):
        ResistiveMHD()
