"""Infinite-space loop quadrature independently checks finite-box vacuum fields."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse.linalg import spsolve
from .electromagnetics import loop_field
from .coupled_benchmarks import aligned_meshes
from .coupling_mesh import CouplingMesh


PROBES = np.asarray(((.75, 0.), (.5, .75), (.75, .5)))


def compact_current(r, z, width=.5):
    """Smooth positive toroidal current [A/m²], supported well inside chamber."""
    return r*np.maximum(1-(r/width)**2, 0)**4*np.maximum(1-(z/width)**2, 0)**4


def infinite_flux(probes=PROBES, order=32):
    """Integrate exact circular-loop Green function against known J dR dz."""
    x, w = np.polynomial.legendre.leggauss(order)
    radius, axial = .25*(x+1), .5*x
    result = np.zeros(len(probes))
    for i, r in enumerate(radius):
        for j, z in enumerate(axial):
            current = compact_current(r, z)*w[i]*w[j]*.125
            result += current*loop_field(probes[:, 0], probes[:, 1], r, z)[0]
    return result


def finite_box_flux(cells, extent=2, probes=PROBES, *, exact_boundary=False):
    """Solve Poisson with zero or independent infinite-space boundary data."""
    mesh, grid = aligned_meshes(cells, extent)
    mapping = CouplingMesh(mesh, grid)
    r, z = grid.r[1:-1, None], grid.z[None, 1:-1]
    rhs = MU0*compact_current(r, z).ravel()*mapping.node_area
    if exact_boundary:
        # Analytic boundary data isolate discretization error from truncating
        # infinite space. Every boundary quadrature point is outside support.
        outer = infinite_flux(np.column_stack((np.full(grid.nz-2, grid.radius), grid.z[1:-1])))
        lower = infinite_flux(np.column_stack((grid.r[1:-1], np.full(grid.nr-2, -grid.half_length))))
        upper = infinite_flux(np.column_stack((grid.r[1:-1], np.full(grid.nr-2, grid.half_length))))
        rows = rhs.reshape(grid.nr-2, grid.nz-2)
        rows[-1] += grid.dz/((grid.radius-grid.dr/2)*grid.dr)*outer
        rows[:, 0] += grid.dr/(grid.r[1:-1]*grid.dz)*lower
        rows[:, -1] += grid.dr/(grid.r[1:-1]*grid.dz)*upper
    solution = spsolve(mapping.k, rhs)
    full = mapping.full(solution)
    ri, zi = probes[:, 0]/grid.dr, (probes[:, 1]+grid.half_length)/grid.dz
    if not np.allclose(ri, np.rint(ri)) or not np.allclose(zi, np.rint(zi)):
        raise ValueError("Verification probes must coincide with grid nodes")
    measured = full[np.rint(ri).astype(int), np.rint(zi).astype(int)]
    return measured, float(np.linalg.norm((mapping.k@solution-rhs))/np.linalg.norm(rhs))
