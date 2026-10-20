"""Variational boundary traction and Poynting flux from the actual CT electric field."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .magnetic_quadrature import stiffness_action


@dataclass
class BoundaryPoynting:
    bz_wall: np.ndarray
    br_bottom: np.ndarray
    br_top: np.ndarray
    bz_inner: np.ndarray | None = None   # inner wall of an annular chamber

    def face_fields(self, mesh):
        w = (mesh.radial_centroid-mesh.rf[:-1])/mesh.dr
        return (.5*(self.bz_wall[:-1]+self.bz_wall[1:]),
                (1-w)*self.br_bottom[:-1]+w*self.br_bottom[1:],
                (1-w)*self.br_top[:-1]+w*self.br_top[1:],
                None if self.bz_inner is None else .5*(self.bz_inner[:-1]+self.bz_inner[1:]))

    def energy_fluxes(self, mesh, electric):
        radial = .5*(electric[-1, :-1]*self.bz_wall[:-1]+electric[-1, 1:]*self.bz_wall[1:])/MU0
        ends = []
        for j, field in ((0, self.br_bottom), (-1, self.br_top)):
            ends.append(-(mesh.rf[:-1]*electric[:-1, j]*field[:-1]+
                          mesh.rf[1:]*electric[1:, j]*field[1:])/(2*mesh.r*MU0))
        inner = None
        if self.bz_inner is not None:   # +r Poynting flux through the inner wall faces
            inner = .5*(electric[0, :-1]*self.bz_inner[:-1]+electric[0, 1:]*self.bz_inner[1:])/MU0
        return radial, ends[0], ends[1], inner


def variational_trace(mapping, full_total, stiffness_current):
    """Boundary load conjugate to psi in vacuum minus fluid magnetic energy.

    Interior rows are deliberately not turned into a boundary correction.
    At a geometric corner the single weak load fixes a weighted sum of two
    traces; a local minimum-norm projection of independent exterior derivatives
    resolves that split consistently for smooth fields.
    """
    mesh = mapping.fluid
    total = mapping.fluid_flux(full_total)
    load = mapping.fluid_flux(mapping.full(stiffness_current))-stiffness_action(mesh, total)
    bz = -load[-1].copy()/mesh.dz
    bottom = -load[:, 0].copy()/mesh.dr
    top = load[:, -1].copy()/mesh.dr
    i0 = mapping.r_offset
    inner = load[0].copy()/mesh.dz if i0 else None   # the inner wall's outward normal is -r
    if not i0:
        bottom[0] = top[0] = 0.
    # The noncorner boundary loads already have full dual boundary lengths.
    # At each wall/end corner both boundary lengths are half an interval; the corner load fixes
    # a*Bz + b*Br and the minimum-norm split moves Bz by -sign(a)*correction, Br by direction*correction.
    corners = [(-1, i0+mesh.nr, -1., bz)]+([(0, i0, 1., inner)] if i0 else [])
    for row, i, wall_sign, wall in corners:
        radius = mesh.rf[row]
        for j, offset, direction, field in ((0, mapping.z_offset, -1, bottom),
                                           (-1, mapping.z_offset+mesh.nz, 1, top)):
            out = 1 if wall_sign < 0 else -1   # one-sided stencil into the vacuum beyond this wall
            bz_guess = out*(-3*full_total[i, offset]+4*full_total[i+out, offset]-full_total[i+2*out, offset])/(2*mesh.dr*radius)
            br_guess = -direction*(-3*full_total[i, offset]+4*full_total[i, offset+direction]-
                                   full_total[i, offset+2*direction])/(2*mesh.dz*radius)
            wall_load, end_load = wall_sign*.5*mesh.dz*bz_guess, direction*.5*mesh.dr*br_guess
            correction = (load[row, j]-wall_load-end_load)/(.5*(mesh.dr+mesh.dz))
            wall[j] = bz_guess+wall_sign*correction
            field[row] = br_guess+direction*correction
    return BoundaryPoynting(bz, bottom, top, inner)
