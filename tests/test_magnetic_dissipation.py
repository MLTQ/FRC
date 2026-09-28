"""Vanishing numerical diffusion, work partition and equivalent implicit solves."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.coupled_benchmarks import aligned_meshes, compact_flux
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.coupled_hydro import gas_from_primitive, gas_primitive
from plasma_bench.coupled_magnetic import MagneticStepper
from plasma_bench.magnetic_dissipation import shock_resistance


def test_smooth_compression_has_h_squared_diffusion_and_no_vacuum_resistance():
    values = []
    for n in (8, 16):
        mesh, grid = aligned_meshes(n)
        mapping = CouplingMesh(mesh, grid, boundary="multipole")
        gas = gas_from_primitive(mesh, 1., 1., vr=-.2*mesh.radial_centroid[:, None], vz=-.1*mesh.z[None, :])
        q = gas_primitive(mesh, gas, 5/3)
        v = np.concatenate((q[..., 1].ravel(), q[..., 2].ravel()))
        resistance, eta = shock_resistance(mapping, q, np.zeros(grid.shape), v, 1., 5/3)
        expected = MU0*.5*mesh.dr**2  # div v = 2*(-.2)+(-.1), including axis cell.
        assert_allclose(eta, expected, rtol=2e-14)
        assert_allclose(resistance[mapping.active]*MU0*mapping.node_area[mapping.active]/mapping.node_radius[mapping.active], expected, rtol=2e-14)
        assert np.all(resistance[mapping.vacuum] == 0)
        assert np.all(shock_resistance(mapping, q, np.zeros(grid.shape), -v, 1., 5/3)[0] == 0)
        assert np.all(shock_resistance(mapping, q, np.zeros(grid.shape), 0*v, 1., 5/3)[0] == 0)
        values.append(eta)
    assert_allclose(values[0]/values[1], 4., rtol=2e-14)


@pytest.mark.parametrize("boundary", ["fixed", "multipole"])
@pytest.mark.parametrize("resistivity", [0., .01*MU0])
@pytest.mark.parametrize("coefficient", [0., 1.])
def test_picard_matches_direct_and_accounts_for_each_heat_source(boundary, resistivity, coefficient):
    mesh, grid = aligned_meshes(8)
    mapping = CouplingMesh(mesh, grid, boundary=boundary)
    induced = compact_flux(grid.r[:, None], grid.z[None, :], amplitude=.1)
    gas = gas_from_primitive(mesh, 1., 1., vr=-.2*mesh.radial_centroid[:, None])
    ext = .5*np.sqrt(MU0)*grid.r[:, None]**2+np.zeros((1, grid.nz))
    updates = [MagneticStepper(mapping, resistivity, method=method, tolerance=1e-13, shock_diffusion=coefficient).step(gas, induced, ext, 1.001*ext, .001)
               for method in ("direct", "picard")]
    assert_allclose(updates[0].gas, updates[1].gas, rtol=2e-11, atol=2e-13)
    assert_allclose(updates[0].induced, updates[1].induced, rtol=2e-10, atol=1e-17)
    for update in updates:
        old, new = mapping.interior(induced), mapping.interior(update.induced)
        delta_field = np.pi/MU0*(new@(mapping.k@new)-old@(mapping.k@old))
        delta_gas = np.sum((update.gas[..., 3]-gas[..., 3])*mesh.volumes)
        assert_allclose(delta_field+delta_gas, update.source_work, atol=3e-14, rtol=3e-11)
        internal = gas_primitive(mesh, update.gas, 5/3)[..., 3]/(5/3-1)
        delta_internal = np.sum((internal-1/(5/3-1))*mesh.volumes)
        assert_allclose(delta_internal, update.ohmic_heat+update.numerical_heat, atol=3e-14, rtol=3e-10)
        assert (update.ohmic_heat > 0) == (resistivity > 0)
        assert (update.numerical_heat > 0) == (coefficient > 0)
        assert np.linalg.norm((mapping.k@new)[mapping.vacuum]) < 1e-11*np.linalg.norm(mapping.k@new)
