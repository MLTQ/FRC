"""Independent conservative overlap, boundary flux and compression checks."""
import numpy as np
import pytest
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_mass import DualMass
from plasma_bench.bulk_entropy import entropy_rate


def test_overlap_volume_integrals_and_adjoint():
    mesh = MHDMesh(.035, .12, 7, 9)
    mapping = DualMass(mesh)
    rng = np.random.default_rng(348)
    primal, nodal = rng.random((7, 9)), rng.random((8, 9))
    np.testing.assert_allclose(mapping.density(np.ones((7, 9))), 1, atol=4e-16)
    np.testing.assert_allclose(np.sum(mesh.volumes*primal), np.sum(mapping.volumes*mapping.density(primal)), rtol=3e-16)
    np.testing.assert_allclose(np.sum(mesh.volumes*primal*mapping.primal_average(nodal)),
                               np.sum(mapping.volumes*mapping.density(primal)*nodal), rtol=3e-16)
    # Q1 interpolation of a linear function has this independently known mean.
    x = .5*mesh.rf**2
    f = np.broadcast_to(3+7*x[:, None], (8, 9))
    np.testing.assert_allclose(mapping.primal_average(f), np.broadcast_to(3+3.5*(x[:-1]+x[1:])[:, None], (7, 9)))


def test_divergence_commutes_locally_with_general_boundary_fluxes():
    mesh = MHDMesh(.035, .12, 7, 9)
    mapping = DualMass(mesh)
    rng = np.random.default_rng(845)
    fr = rng.normal(size=(8, 9))
    fz = rng.normal(size=(7, 9))
    fz = np.column_stack((fz, fz[:, 0]))
    primal = -(np.diff(fr, axis=0)+np.diff(fz, axis=1))/mesh.volumes
    dr, dz = mapping.fluxes(fr, fz)
    actual = -mapping.geometry.divergence(dr, dz)
    np.testing.assert_allclose(actual, mapping.density(primal), atol=4e-9, rtol=3e-14)
    assert abs(np.sum(mapping.volumes*actual)+np.sum(fr[-1]-fr[0])) < 5e-15
    with pytest.raises(ValueError, match='periodic'):
        mapping.fluxes(fr, np.ones((7, 10))*np.arange(10))


def test_affine_contravariant_flux_matches_exact_dual_face_integrals():
    mesh = MHDMesh(.04, .15, 9, 11)
    mapping = DualMass(mesh)
    x = .5*mesh.rf**2
    # div(F) = d(Fx)/dx = 2.3; volume measure is 2pi dx dz.
    fr = np.broadcast_to(2*np.pi*mesh.dz*(2.3*x[:, None]+.7), (10, 11))
    fz = np.broadcast_to(2*np.pi*np.diff(x)[:, None]*.4, (9, 12))
    dr, dz = mapping.fluxes(fr, fz)
    np.testing.assert_allclose(dr, np.broadcast_to(2*np.pi*mesh.dz*(2.3*mapping.geometry.dual_x_faces[:, None]+.7), dr.shape), rtol=4e-16)
    np.testing.assert_allclose(dz, np.broadcast_to(2*np.pi*np.diff(mapping.geometry.dual_x_faces)[:, None]*.4, dz.shape), rtol=8e-16)
    np.testing.assert_allclose(mapping.geometry.divergence(dr, dz), 2.3, rtol=3e-10)


@pytest.mark.parametrize('order', [1, 2])
def test_constant_specific_entropy_compresses_with_shared_mass(order):
    mesh = MHDMesh(.04, .15, 9, 11)
    mapping = DualMass(mesh)
    rng = np.random.default_rng(921)
    n = mapping.density(1+rng.random((9, 11)))
    fr = rng.normal(size=(10, 11)); fr[[0, -1]] = 0
    fz = rng.normal(size=(9, 11)); fz = np.column_stack((fz, fz[:, 0]))
    dr, dz = mapping.fluxes(fr, fz)
    dn = -mapping.geometry.divergence(dr, dz)
    du = entropy_rate(mapping, n, .031*n, dr, dz, order=order)
    np.testing.assert_allclose(du, .031*dn, atol=3e-10, rtol=4e-14)
    assert abs(np.sum(mapping.volumes*du)) < 2e-16
