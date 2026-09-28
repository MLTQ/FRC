"""Asymmetric smooth current and independent infinite-space verification."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse.linalg import spsolve
from .coupled_benchmarks import aligned_meshes
from .coupling_mesh import CouplingMesh
from .electromagnetics import loop_field


@dataclass(frozen=True)
class CurrentPacket:
    radius_center: float = .5
    axial_center: float = .35
    radial_width: float = .25
    axial_width: float = .4

    def value(self, r, z):
        return r*np.maximum(1-((r-self.radius_center)/self.radial_width)**2, 0)**4*\
            np.maximum(1-((z-self.axial_center)/self.axial_width)**2, 0)**4

    def reference(self, probes, order=32):
        """Quadrature of exact circular-loop (psi, Br, Bz), without discretized K."""
        probes = np.asarray(probes)
        x, w = np.polynomial.legendre.leggauss(order)
        radial = self.radius_center+self.radial_width*x
        axial = self.axial_center+self.axial_width*x
        result = np.zeros((3, len(probes)))
        for i, r in enumerate(radial):
            for j, z in enumerate(axial):
                current = self.value(r, z)*w[i]*w[j]*self.radial_width*self.axial_width
                result += current*np.asarray(loop_field(probes[:, 0], probes[:, 1], r, z))
        return result


INTERIOR_PROBES = np.asarray(((.875, .25), (.5, -.5), (.25, .875)))
EXTERIOR_PROBES = np.asarray(((.6, 2.5), (2.5, 0.), (2.5, .8), (.4, -3.), (0., 3.)))


def static_packet(cells, *, extent=2, modes=12, packet=None):
    """Solve the current packet with infinite magnetic exterior, no fluid update."""
    packet = packet or CurrentPacket()
    mesh, grid = aligned_meshes(cells, extent)
    mapping = CouplingMesh(mesh, grid, boundary="multipole", multipole_modes=modes)
    current = np.zeros(mapping.node_count)
    current[:mapping.grid_node_count] = packet.value(grid.r[1:-1, None], grid.z[None, 1:-1]).ravel()
    rhs = MU0*current*mapping.node_area
    unknowns = spsolve(mapping.k, rhs)
    return mapping, unknowns, rhs


def nodal_probes(mapping, unknowns, probes=INTERIOR_PROBES):
    """Sample actual grid nodes, with no interpolation error in the comparison."""
    grid = mapping.field
    ri, zi = probes[:, 0]/grid.dr, (probes[:, 1]+grid.half_length)/grid.dz
    if not np.allclose(ri, np.rint(ri)) or not np.allclose(zi, np.rint(zi)):
        raise ValueError("Probe positions must coincide with actual field nodes")
    return mapping.full(unknowns)[np.rint(ri).astype(int), np.rint(zi).astype(int)]
