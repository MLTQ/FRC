import numpy as np
from numpy.testing import assert_allclose
from plasma_bench.geometry import Grid
from plasma_bench.topology import analyze_flux
from plasma_bench.diagnostics import longest_sampled_interval


def test_uniform_and_uniform_reversed_field_have_no_closed_flux():
    g = Grid(.12, .18, 61, 101)
    r, z = g.mesh
    for sign in (-1, 1):
        d = analyze_flux(g, sign*r*r)
        assert not d["closed_flux"] and not d["o_points"]


def test_manufactured_frc_flux_axis_nulls_and_volume():
    g = Grid(.12, .18, 101, 151)
    r, z = g.mesh
    a, b = .08, .12
    psi = r*r*(1-r*r/a**2-z*z/b**2)
    d = analyze_flux(g, psi)
    assert d["closed_flux"] and len(d["o_points"]) == 1
    assert_allclose(d["o_points"][0]["r_m"], a/np.sqrt(2), rtol=.002)
    assert_allclose([x["z_m"] for x in d["axis_nulls"]], [-b, b], atol=g.dz)
    assert_allclose(d["separatrix_radius_m"], a, rtol=.01)
    assert_allclose(d["separatrix_volume_m3"], 4*np.pi*a*a*b/3, rtol=.015)


def test_off_axis_saddle_classification():
    g = Grid(.12, .18, 81, 121)
    r, z = g.mesh
    # Local saddle is interior; this fixture only tests classification, not a physical boundary.
    psi = (r-.06)**2-z*z
    d = analyze_flux(g, psi)
    assert len(d["x_points"]) == 1 and not d["closed_flux"]


def test_lifetime_is_sampled_lower_bound():
    assert longest_sampled_interval(np.arange(5.), [False, True, True, False, True]) == 1
    assert longest_sampled_interval(np.arange(3.), [False, True, False]) == 0
