import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.constants import mu_0
from plasma_bench.electromagnetics import loop_field, VacuumField
from plasma_bench.geometry import Grid
from plasma_bench.coils import Coil
from plasma_bench.waveforms import DC, Gaussian
from plasma_bench.numerics import field_from_flux, cell_divergence


def test_axis_exact_and_finite():
    a = 0.07
    z = np.linspace(-0.4, 0.4, 51)
    psi, br, bz = loop_field(0, z, a)
    assert_allclose(bz, mu_0*a*a/(2*(a*a+z*z)**1.5), rtol=3e-15)
    assert np.all(psi == 0) and np.all(br == 0)


@pytest.mark.parametrize("r,z", [(0.02, 0.03), (0.12, -0.04), (1e-8, 0.06), (0.4, 1.0)])
def test_against_independent_biot_savart_quadrature(r, z):
    a = 0.07
    def distance(phi):
        return (r*r+a*a-2*r*a*np.cos(phi)+z*z)**1.5
    br_ref = mu_0/(4*np.pi)*quad(lambda p: a*z*np.cos(p)/distance(p), 0, 2*np.pi, epsabs=1e-12)[0]
    bz_ref = mu_0/(4*np.pi)*quad(lambda p: (a*a-a*r*np.cos(p))/distance(p), 0, 2*np.pi, epsabs=1e-12)[0]
    _, br, bz = loop_field(r, z, a)
    assert_allclose([br, bz], [br_ref, bz_ref], rtol=2e-8, atol=2e-19)


def test_flux_equals_integrated_axial_field():
    for r in (1e-7, 0.015, 0.05, 0.12):
        p, _, _ = loop_field(r, 0.03, 0.07)
        integrated = quad(lambda x: x*float(loop_field(x, 0.03, 0.07)[2]), 0, r, epsabs=1e-22)[0]
        assert_allclose(p, integrated, rtol=3e-10, atol=1e-28)


def test_symmetry_superposition_and_far_dipole():
    p, br, bz = loop_field(0.03, 0.04, 0.07)
    p2, br2, bz2 = loop_field(0.03, -0.04, 0.07)
    assert_allclose([p, br, bz], [p2, -br2, bz2])
    _, _, far = loop_field(0, 7.0, 0.07)
    assert_allclose(far, mu_0*0.07**2/(2*7.0**3), rtol=2e-4)
    g = Grid(0.04, 0.1, 21, 41)
    c = Coil(0.07, 0, DC(10), turns=3)
    v = VacuumField(g, (c, c), 0.01).at(0)
    unit = loop_field(*g.mesh, 0.07)
    assert_allclose(v.bz, 0.01+60*unit[2])
    assert np.all(v.e_theta == 0)


def test_faraday_and_polarity():
    g = Grid(0.04, 0.1, 21, 41)
    w = Gaussian(100, 2e-5, 1e-5)
    v = VacuumField(g, (Coil(0.07, 0, w, polarity=-1),))
    t, dt = 1e-5, 1e-10
    f = v.at(t)
    derivative = (v.at(t+dt).psi-v.at(t-dt).psi)/(2*dt)
    assert_allclose(f.dpsi_dt, derivative, rtol=2e-10, atol=1e-16)
    assert_allclose(g.r[:, None]*f.e_theta, -f.dpsi_dt)


def test_flux_reconstruction_converges_and_discrete_divergence():
    errors = []
    for n in (21, 41, 81):
        grid = Grid(0.04, 0.1, n, 2*n-1)
        p, br, bz = loop_field(*grid.mesh, 0.07)
        brh, bzh = field_from_flux(grid, p)
        errors.append(np.max(np.abs(bzh[2:-2, 2:-2]-bz[2:-2, 2:-2])))
        scale = np.max(np.abs(bz))/min(grid.dr, grid.dz)
        assert np.max(np.abs(cell_divergence(grid, p))) < 1e-12*scale
        assert_allclose(bzh[0], bz[0], rtol=1e-5 if n == 21 else 1e-6)
    assert errors[0]/errors[1] > 3.4
    assert errors[1]/errors[2] > 3.4


def test_uniform_flux_axis_and_singular_wire():
    g = Grid()
    r, z = g.mesh
    br, bz = field_from_flux(g, 0.03*r*r/2)
    assert_allclose(bz, 0.03, rtol=2e-13)
    assert_allclose(br, 0, atol=1e-14)
    with pytest.raises(ValueError, match="singular"):
        loop_field(0.07, 0, 0.07)
