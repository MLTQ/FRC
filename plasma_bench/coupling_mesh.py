"""Aligned fluid cells, current-bearing nodes and an exact elliptic exterior."""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from .numerics import stiffness, expand
from .open_boundary import MultipoleBoundary


def _linear_weights(centers, location):
    high = int(np.clip(np.searchsorted(centers, location), 1, len(centers)-1))
    low = high-1
    upper = (location-centers[low])/(centers[high]-centers[low])
    return ((low, 1-upper), (high, upper))


class CouplingMesh:
    """Fluid exists only in its declared cylinder; no exterior fluid arrays.

    The aligned field grid extends beyond both radius and ends. The induced
    exterior equation is K u=0, with either a fixed far boundary or an infinite
    multipole exterior. The fluid interface is fixed, not a moving free surface.
    """
    def __init__(self, fluid, field, *, boundary="fixed", multipole_modes=12):
        self.fluid, self.field = fluid, field
        if not np.isclose(fluid.dr, field.dr, rtol=1e-12, atol=0) or not np.isclose(fluid.dz, field.dz, rtol=1e-12, atol=0):
            raise ValueError("Fluid faces and field nodes require matching dr and dz")
        offset = (field.half_length-fluid.half_length)/field.dz
        if not np.isclose(offset, round(offset), rtol=0, atol=1e-10):
            raise ValueError("Fluid end faces must align with field nodes")
        self.z_offset = int(round(offset))
        if field.radius < fluid.radius+field.dr*(2-1e-12) or self.z_offset < 2:
            raise ValueError("At least two exterior field intervals are required radially and beyond both ends")
        inner = fluid.inner_radius/field.dr
        if not np.isclose(inner, round(inner), rtol=0, atol=1e-10):
            raise ValueError("An inner chamber wall must lie on a field node row")
        self.r_offset = int(round(inner))      # field row of the first fluid face (0: the axis)
        if 0 < self.r_offset < 2:
            raise ValueError("At least two vacuum field intervals are required inside an inner chamber wall")
        if boundary not in ("fixed", "multipole"):
            raise ValueError("Magnetic boundary must be 'fixed' or 'multipole'")
        if boundary == "multipole" and min(field.radius, field.half_length) <= np.hypot(fluid.radius, fluid.half_length):
            raise ValueError("Multipole matching boundary must enclose the entire fluid source sphere")
        self.open_boundary = MultipoleBoundary(field, multipole_modes) if boundary == "multipole" else None
        self.k = self.open_boundary.k if self.open_boundary else stiffness(field)
        self.grid_node_count = (field.nr-2)*(field.nz-2)
        self.node_count, self.cell_count = self.k.shape[0], fluid.nr*fluid.nz
        self.node_radius = np.pad(np.broadcast_to(field.r[1:-1, None], (field.nr-2, field.nz-2)).ravel(),
                                  (0, self.node_count-self.grid_node_count))
        self.node_area = np.zeros(self.node_count)
        rr, rc, rv, zr, zc, zv, hr, hc, hv = ([] for _ in range(9))
        for i in range(0 if self.r_offset else 1, fluid.nr+1):   # the axis row is fixed, an inner wall is not
            for j in range(fluid.nz+1):
                node = (self.r_offset+i-1)*(field.nz-2)+(self.z_offset+j-1)
                radial_cells = [k for k in (i-1, i) if 0 <= k < fluid.nr]
                axial_cells = [k for k in (j-1, j) if 0 <= k < fluid.nz]
                cells = [a*fluid.nz+b for a in radial_cells for b in axial_cells]
                self.node_area[node] = len(cells)*fluid.dr*fluid.dz/4
                for cell in cells:
                    hr.append(cell); hc.append(node); hv.append(1/len(cells))
                rw = _linear_weights(fluid.radial_centroid, fluid.rf[i])
                zw = _linear_weights(fluid.z, fluid.zf[j])
                for a, wa in rw:
                    for b, wb in zw:
                        cell = a*fluid.nz+b
                        if 0 < i < fluid.nr:  # Impermeable radial chamber walls.
                            rr.append(node); rc.append(cell); rv.append(wa*wb)
                        if 0 < j < fluid.nz:  # Impermeable chamber ends.
                            zr.append(node); zc.append(cell); zv.append(wa*wb)
        shape = (self.node_count, self.cell_count)
        self.interpolate_r = coo_matrix((rv, (rr, rc)), shape=shape).tocsr()
        self.interpolate_z = coo_matrix((zv, (zr, zc)), shape=shape).tocsr()
        self.heat_to_cells = coo_matrix((hv, (hr, hc)), shape=shape[::-1]).tocsr()
        self.active = self.node_area > 0
        self.vacuum = ~self.active
        self.active_diagonal = diags(self.active.astype(float), format="csc")
        self.vacuum_rows = diags(self.vacuum.astype(float), format="csc")@self.k
        self._vacuum_factor = splu(self.k[self.vacuum][:, self.vacuum].tocsc())
        self._vacuum_cross = self.k[self.vacuum][:, self.active]

    def interior(self, flux):
        flux = np.asarray(flux, dtype=float)
        if flux.shape != self.field.shape or not np.all(np.isfinite(flux)):
            raise ValueError("Flux must be finite on the full field-node grid")
        return self.open_boundary.pack(flux) if self.open_boundary else flux[1:-1, 1:-1].ravel().copy()

    def grid_vector(self, flux):
        """Physical grid values, padded by zero for non-nodal modal unknowns."""
        flux = np.asarray(flux, dtype=float)
        if flux.shape != self.field.shape or not np.all(np.isfinite(flux)):
            raise ValueError("Flux must be finite on the full field-node grid")
        return np.pad(flux[1:-1, 1:-1].ravel(), (0, self.node_count-self.grid_node_count))

    def full(self, interior):
        return self.open_boundary.expand(interior) if self.open_boundary else expand(self.field, interior)

    def fluid_flux(self, full_flux):
        self.grid_vector(full_flux)
        return np.asarray(full_flux)[self.r_offset:self.r_offset+self.fluid.nr+1, self.z_offset:self.z_offset+self.fluid.nz+1]

    def harmonic_extension(self, flux):
        """Keep the prescribed plasma-node flux and solve only exterior vacuum."""
        return self.full(self.complete_vacuum(self.grid_vector(flux)))

    def complete_vacuum(self, unknowns):
        """Retain plasma unknowns and solve cached exterior constraints."""
        u = np.asarray(unknowns, float).copy()
        if u.shape != (self.node_count,) or not np.all(np.isfinite(u)):
            raise ValueError("Finite magnetic unknowns must match the coupled mesh")
        u[self.vacuum] = self._vacuum_factor.solve(-self._vacuum_cross@u[self.active])
        return u
