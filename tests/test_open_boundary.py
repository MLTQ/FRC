import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.constants import mu_0 as MU0
from plasma_bench.geometry import Grid
from plasma_bench.open_boundary import multipoles, exterior_energy, MultipoleBoundary
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.coupled_benchmarks import aligned_meshes, coil_pulse
from plasma_bench.open_benchmarks import CurrentPacket, static_packet, nodal_probes, INTERIOR_PROBES, EXTERIOR_PROBES


def test_dipole_and_derivatives_match_closed_form_including_axis():
    r, z, scale = np.asarray([.4, 1., 2., 0.]), np.asarray([3., -2., .5, -3.]), 2.
    radius = np.hypot(r, z); amplitude = scale/np.sqrt(2)
    psi, dr, dz = multipoles(r, z, scale, 1)
    assert_allclose(psi[:, 0], amplitude*r*r/radius**3, rtol=3e-14, atol=1e-15)
    assert_allclose(dr[:, 0], amplitude*(2*r/radius**3-3*r**3/radius**5), rtol=3e-14, atol=1e-15)
    assert_allclose(dz[:, 0], -3*amplitude*r*r*z/radius**5, rtol=3e-14, atol=1e-15)
    boundary = MultipoleBoundary(Grid(2, 2, 9, 17), 1)
    flux, br, bz = boundary.outside_field([1.], r, z)
    assert_allclose(bz[-1], 2*amplitude/abs(z[-1])**3, rtol=1e-14)
    assert flux[-1] == br[-1] == 0


def test_exterior_energy_against_independent_dipole_volume_integral():
    grid = Grid(2, 3, 9, 17)
    # Integrate analytic dipole |grad psi|²/r in polar coordinates, performing
    # the entire R integral analytically from the rectangular boundary to infinity.
    amplitude = min(grid.radius, grid.half_length)/np.sqrt(2)
    angle = np.arctan(grid.radius/grid.half_length)
    def integrand(theta):
        s, c = np.sin(theta), abs(np.cos(theta))
        start = min(grid.radius/s, grid.half_length/c) if c > 0 else grid.radius/s
        return amplitude**2*(s**3+4*s*c*c)/(3*start**3)
    exact = sum(quad(integrand, a, b, epsabs=1e-13)[0] for a, b in ((0, angle), (angle, np.pi-angle), (np.pi-angle, np.pi)))
    measured = exterior_energy(grid, 12)
    assert_allclose(measured[0, 0], exact, rtol=2e-13)
    assert_allclose(exterior_energy(grid, 12, quadrature_order=128), measured, rtol=1e-11, atol=1e-13)
    assert np.linalg.eigvalsh(measured).min() > 0


def test_modal_flux_is_harmonic_and_decays_to_infinity():
    r, z = np.asarray([1., 2., 3.]), np.asarray([3., 2., -1.])
    h = 2e-4
    psi, dr, _ = multipoles(r, z, 2., 8)
    radial_lap = (multipoles(r+h, z, 2., 8)[0]-2*psi+multipoles(r-h, z, 2., 8)[0])/h**2
    axial_lap = (multipoles(r, z+h, 2., 8)[0]-2*psi+multipoles(r, z-h, 2., 8)[0])/h**2
    assert np.max(abs(radial_lap+axial_lap-dr/r[:, None])) < 1e-7
    assert_allclose(multipoles(2*r, 2*z, 2., 8)[0], psi/2.**np.arange(1, 9), rtol=2e-14)


def test_asymmetric_current_matches_independent_infinite_space_field():
    packet = CurrentPacket()
    exact = packet.reference(INTERIOR_PROBES)[0]
    assert_allclose(packet.reference(EXTERIOR_PROBES, 48), packet.reference(EXTERIOR_PROBES), rtol=2e-12, atol=1e-23)
    errors = []
    for cells in (16, 32):
        mapping, unknowns, rhs = static_packet(cells)
        errors.append(np.linalg.norm(nodal_probes(mapping, unknowns)-exact)/np.linalg.norm(exact))
        assert np.linalg.norm(mapping.k@unknowns-rhs)/np.linalg.norm(rhs) < 1e-11
        assert (mapping.k-mapping.k.T).nnz == 0
        recovered = mapping.interior(mapping.full(unknowns))
        assert np.linalg.norm(recovered-unknowns) < 5e-13*np.linalg.norm(unknowns)
    assert errors[1] < errors[0]/3.5


def test_retaining_more_multipoles_resolves_asymmetric_exterior_field():
    exact = CurrentPacket().reference(EXTERIOR_PROBES)
    errors = []
    for modes in (1, 2, 4, 12):
        mapping, unknowns, _ = static_packet(32, modes=modes)
        field = np.asarray(mapping.open_boundary.outside_field(unknowns[mapping.grid_node_count:], EXTERIOR_PROBES[:, 0], EXTERIOR_PROBES[:, 1]))
        errors.append(np.linalg.norm(field[1:]-exact[1:])/np.linalg.norm(exact[1:]))
    assert errors[1] < errors[0]/2
    assert errors[2] < errors[1]/3
    assert errors[-1] < .003


def test_open_boundary_pulse_is_insensitive_to_matching_surface_and_conserves_work():
    rows = []
    for extent in (2, 3):
        solver = coil_pulse(8, extent=extent, vacuum_boundary="multipole")
        solver.advance(.1)
        d = solver.diagnostics(); rows.append(d)
        assert d["beyond_grid_induced_energy_j"] > 0
        assert abs(d["response_energy_balance_relative"]) < 2e-12
        assert d["open_boundary_relative_stationarity_residual"] < 1e-10
        assert d["max_vacuum_induced_current_residual_a_m2"] < 1e-7
    assert abs(rows[0]["source_work_j"]/rows[1]["source_work_j"]-1) < 2e-4
    assert abs(rows[0]["maximum_speed_m_s"]/rows[1]["maximum_speed_m_s"]-1) < 2e-4


def test_open_boundary_input_contracts():
    mesh, grid = aligned_meshes(8)
    with pytest.raises(ValueError, match="source sphere"):
        CouplingMesh(mesh, Grid(1.25, 1.25, 11, 21), boundary="multipole")
    with pytest.raises(ValueError, match="integer"):
        MultipoleBoundary(grid, 0)
    boundary = MultipoleBoundary(grid, 4)
    arbitrary = np.zeros(grid.shape); arbitrary[-1, 3] = 1
    with pytest.raises(ValueError, match="representable"):
        boundary.pack(arbitrary)
    with pytest.raises(ValueError, match="inside"):
        boundary.outside_field(np.zeros(4), .5, .5)
