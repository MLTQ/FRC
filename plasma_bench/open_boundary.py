"""Variational decaying multipoles for the infinite axisymmetric magnetic exterior."""
import numpy as np
from scipy.special import eval_legendre
from scipy.sparse import bmat, coo_matrix, csc_matrix
from .numerics import stiffness


def multipoles(r, z, scale, modes):
    """Return psi, dpsi/dr, dpsi/dz for regular-axis decaying vacuum modes.

    Column l-1 is (scale/R)^l (1-x²) P_l'(x)/sqrt(l(l+1)), l>=1,
    where R=hypot(r,z), x=z/R. Evaluate only outside the source-enclosing sphere.
    """
    r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
    radius = np.hypot(r, z)
    if np.any(r < 0) or np.any(radius <= 0) or not np.all(np.isfinite(radius)):
        raise ValueError("Multipole coordinates require r>=0 and finite R>0")
    x, sine = z/radius, r/radius
    flux, radial, axial = [], [], []
    for ell in range(1, modes+1):
        p = eval_legendre(ell, x)
        f = ell*(eval_legendre(ell-1, x)-x*p)
        derivative = -ell*(ell+1)*p
        factor = (scale/radius)**ell/np.sqrt(ell*(ell+1))
        flux.append(factor*f)
        radial.append(-factor/radius*sine*(ell*f+x*derivative))
        axial.append(factor/radius*(-ell*x*f+(1-x*x)*derivative))
    return tuple(np.stack(values, axis=-1) for values in (flux, radial, axial))


def exterior_energy(grid, modes, *, quadrature_order=None):
    """Exact-to-quadrature exterior Dirichlet energy via the inner boundary."""
    order = quadrature_order or max(64, 4*modes+16)
    x, w = np.polynomial.legendre.leggauss(order)
    scale = min(grid.radius, grid.half_length)
    matrix = np.zeros((modes, modes))
    # Normals point out of the finite box, opposite to the exterior's normal.
    edges = ((np.full(order, grid.radius), grid.half_length*x, grid.half_length*w, 0, 1),
             (.5*grid.radius*(x+1), np.full(order, grid.half_length), .5*grid.radius*w, 1, 1),
             (.5*grid.radius*(x+1), np.full(order, -grid.half_length), .5*grid.radius*w, 1, -1))
    for r, z, weight, component, sign in edges:
        psi, dr, dz = multipoles(r, z, scale, modes)
        normal = sign*(dr if component == 0 else dz)
        matrix -= psi.T@((weight/r)[:, None]*normal)
    skew = np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)
    if skew > 1e-10:
        raise ValueError("Exterior energy quadrature is not sufficiently reciprocal")
    matrix = .5*(matrix+matrix.T)
    if np.linalg.eigvalsh(matrix)[0] <= 0:
        raise ValueError("Exterior multipole energy must be positive definite")
    return matrix


class MultipoleBoundary:
    """Unknowns are interior nodal flux followed by exterior multipole amplitudes.

    A symmetric energy matrix combines the interior discrete field and the
    infinite exterior. Boundary nodes share the multipole trace; no zero-flux
    artificial wall or exterior current is imposed at this matching surface.
    """
    def __init__(self, grid, modes=12):
        if not isinstance(modes, int) or isinstance(modes, bool) or modes < 1:
            raise ValueError("At least one integer multipole mode is required")
        self.grid, self.modes = grid, modes
        self.scale = min(grid.radius, grid.half_length)
        self.interior_count = (grid.nr-2)*(grid.nz-2)
        # Every outer boundary node except the symmetry axis, with no duplicates.
        self.bi = np.concatenate((np.full(grid.nz, grid.nr-1), np.arange(1, grid.nr-1), np.arange(1, grid.nr-1)))
        self.bj = np.concatenate((np.arange(grid.nz), np.zeros(grid.nr-2, int), np.full(grid.nr-2, grid.nz-1)))
        self.trace = multipoles(grid.r[self.bi], grid.z[self.bj], self.scale, modes)[0]
        singular = np.linalg.svd(self.trace, compute_uv=False)
        if len(singular) < modes or singular[-1] < 1e-10*singular[0]:
            raise ValueError("Too many multipoles for the matching boundary grid")
        self._inverse_trace = np.linalg.pinv(self.trace, rcond=1e-12)
        self.exterior = exterior_energy(grid, modes)
        cross, boundary = self._boundary_blocks()
        coupling = cross@self.trace
        lower = self.trace.T@(boundary@self.trace)+self.exterior
        lower = .5*(lower+lower.T)
        self.k = bmat([[stiffness(grid), csc_matrix(coupling)],
                       [csc_matrix(coupling.T), csc_matrix(lower)]], format="csc")

    def _boundary_blocks(self):
        g = self.grid
        boundary_count = len(self.bi)
        ids = np.full(g.shape, -1, int)
        ids[1:-1, 1:-1] = np.arange(self.interior_count).reshape(g.nr-2, g.nz-2)
        ids[self.bi, self.bj] = self.interior_count+np.arange(boundary_count)
        br, bc, bv, cr, cc, cv = ([] for _ in range(6))

        def add_edge(a, b, weight):
            # Only boundary blocks are assembled; K_II remains the existing
            # verified stencil. Axis endpoints have fixed zero flux (id=-1).
            for row, col, value in ((a, a, weight), (b, b, weight), (a, b, -weight), (b, a, -weight)):
                if row < 0 or col < 0:
                    continue
                if row >= self.interior_count and col >= self.interior_count:
                    br.append(row-self.interior_count); bc.append(col-self.interior_count); bv.append(value)
                elif 0 <= row < self.interior_count <= col:
                    cr.append(row); cc.append(col-self.interior_count); cv.append(value)

        for i in range(g.nr-1):
            for j in range(g.nz):
                if j in (0, g.nz-1) or i == g.nr-2:
                    dz = g.dz*(.5 if j in (0, g.nz-1) else 1.)
                    add_edge(ids[i, j], ids[i+1, j], dz/((g.r[i]+g.dr/2)*g.dr))
        for i in range(1, g.nr):
            for j in range(g.nz-1):
                if i == g.nr-1 or j in (0, g.nz-2):
                    dr = g.dr*(.5 if i == g.nr-1 else 1.)
                    add_edge(ids[i, j], ids[i, j+1], dr/(g.r[i]*g.dz))
        return (coo_matrix((cv, (cr, cc)), shape=(self.interior_count, boundary_count)).tocsr(),
                coo_matrix((bv, (br, bc)), shape=(boundary_count, boundary_count)).tocsr())

    def pack(self, flux):
        boundary = flux[self.bi, self.bj]
        coefficients = self._inverse_trace@boundary
        residual = np.linalg.norm(self.trace@coefficients-boundary)
        if residual > 1e-10*max(np.linalg.norm(flux), 1e-30):
            raise ValueError("Induced boundary flux is not representable by the selected multipoles")
        return np.concatenate((flux[1:-1, 1:-1].ravel(), coefficients))

    def expand(self, unknowns):
        result = np.zeros(self.grid.shape)
        result[1:-1, 1:-1] = unknowns[:self.interior_count].reshape(self.grid.nr-2, self.grid.nz-2)
        result[self.bi, self.bj] = self.trace@unknowns[self.interior_count:]
        return result

    def outside_field(self, amplitudes, r, z):
        """Evaluate (psi, Br, Bz) outside/on the matching rectangle, including axis."""
        amplitudes = np.asarray(amplitudes, float)
        r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
        if amplitudes.shape != (self.modes,) or not np.all(np.isfinite(amplitudes)):
            raise ValueError("Finite amplitudes must match the retained multipole count")
        if np.any((r < self.grid.radius)&(abs(z) < self.grid.half_length)):
            raise ValueError("Exterior multipoles cannot evaluate points inside the matching rectangle")
        psi, dr, dz = multipoles(r, z, self.scale, self.modes)
        flux, pr, pz = psi@amplitudes, dr@amplitudes, dz@amplitudes
        br = np.divide(-pz, r, out=np.zeros_like(pz), where=r > 0)
        bz = np.divide(pr, r, out=np.zeros_like(pr), where=r > 0)
        ell = np.arange(1, self.modes+1)
        radius = np.hypot(r, z)
        axis_limit = ((self.scale/radius[..., None])**ell/radius[..., None]**2*
                      np.sqrt(ell*(ell+1))*np.sign(z[..., None])**(ell-1))@amplitudes
        return flux, br, np.where(r == 0, axis_limit, bz)
