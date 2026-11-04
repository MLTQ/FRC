"""Three-dimensional vacuum field of an open wire cage in a grounded sphere."""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import epsilon_0
from scipy.linalg import lstsq


def directions(n, phase=0.):
    """Approximately uniform sphere directions, for independent sampling."""
    k = np.arange(n)+.5
    z = 1-2*k/n
    angle = k*np.pi*(3-np.sqrt(5))+phase
    return np.column_stack((np.sqrt(1-z*z)*np.cos(angle),
                            np.sqrt(1-z*z)*np.sin(angle), z))


@dataclass
class Cage:
    vertices: np.ndarray
    edges: np.ndarray
    wire_radius: float
    wall_radius: float
    cap_directions: np.ndarray

    def wire_distance(self, points):
        """Distance to the union of finite edge axes (rounded wire junctions)."""
        points = np.asarray(points, float).reshape(-1, 3)
        a, b = self.vertices[self.edges[:, 0]], self.vertices[self.edges[:, 1]]
        u = b-a
        result = np.empty(len(points))
        for start in range(0, len(points), 256):
            x = points[start:start+256, None, :]-a
            t = np.clip(np.sum(x*u, axis=-1)/np.sum(u*u, axis=-1), 0., 1.)
            result[start:start+256] = np.sqrt(np.min(np.sum((x-t[..., None]*u)**2, axis=-1), axis=1))
        return result


def truncated_cage(diameter=.15, wire_diameter=.002, wall_radius=.2):
    """Truncate all directed edges of an icosahedron; retain all 90 wires."""
    if not 0 < wire_diameter < diameter < 2*wall_radius:
        raise ValueError('Require 0 < wire diameter < cage diameter < chamber diameter')
    golden = (1+np.sqrt(5))/2
    ico = np.array([(0, a, b*golden) for a in (-1, 1) for b in (-1, 1)]
                   + [(b*golden, 0, a) for a in (-1, 1) for b in (-1, 1)]
                   + [(a, b*golden, 0) for a in (-1, 1) for b in (-1, 1)], float)
    dist = np.linalg.norm(ico[:, None]-ico[None, :], axis=-1)
    directed = np.argwhere(np.isclose(dist, np.min(dist[dist > 0])))
    vertices = np.array([(2*ico[i]+ico[j])/3 for i, j in directed])
    vertices *= (diameter-wire_diameter)/(2*np.linalg.norm(vertices[0]))
    dist = np.linalg.norm(vertices[:, None]-vertices[None, :], axis=-1)
    edges = np.argwhere(np.triu(np.isclose(dist, np.min(dist[dist > 0])), 1))
    if len(vertices) != 60 or len(edges) != 90 or not np.all(np.bincount(edges.ravel()) == 3):
        raise ValueError('Invalid truncated-icosahedron connectivity')
    return Cage(vertices, edges, wire_diameter/2, wall_radius,
                ico/np.linalg.norm(ico, axis=1)[:, None])


@dataclass
class Sources:
    """Harmonic basis sources entirely inside metal, not physical plasma."""
    starts: np.ndarray
    ends: np.ndarray
    points: np.ndarray
    wall_radius: float
    quadrature: int = 8

    @property
    def count(self):
        return len(self.starts)+len(self.points)

    def images(self, x, y, field):
        norm = np.linalg.norm(y, axis=1)
        live = norm > 0
        potential = np.full((len(x), len(y)), -1/self.wall_radius)
        electric = np.zeros((len(x), len(y), 3)) if field else None
        image = self.wall_radius**2*y[live]/norm[live, None]**2
        factor = self.wall_radius/norm[live]
        dx = x[:, None, :]-image
        radius = np.linalg.norm(dx, axis=-1)
        potential[:, live] = -factor/radius
        if field:
            electric[:, live] = -factor[None, :, None]*dx/radius[..., None]**3
        return potential, electric

    def kernel(self, x, field=False):
        """Potential per strength q/(4*pi*eps0), optionally its negative gradient.

        Direct finite-line integrals are analytic. Sphere image terms use
        Gauss quadrature; no softened Coulomb kernel is used in vacuum.
        """
        x = np.asarray(x, float).reshape(-1, 3)
        count = len(self.starts)
        potential = np.zeros((len(x), self.count))
        electric = np.zeros((len(x), self.count, 3)) if field else None
        if count:
            delta = self.ends-self.starts
            length = np.linalg.norm(delta, axis=1)
            u = delta/length[:, None]
            dx = x[:, None, :]-self.starts
            t = np.sum(dx*u, axis=-1)
            perp = dx-t[..., None]*u
            r2 = np.sum(perp*perp, axis=-1)
            r = np.sqrt(r2)
            on_axis = r < 1e-12*length
            singular = on_axis & (t >= 0) & (t <= length)
            if singular.any():
                raise ValueError('Evaluation point lies on an interior source line')
            safe_r = np.where(on_axis, 1., r)
            value = (np.arcsinh((length-t)/safe_r)+np.arcsinh(t/safe_r))/length
            near = np.where(t < 0, -t, t-length)
            value[on_axis] = np.broadcast_to(np.log1p(length/np.maximum(near, 1e-300))/length, t.shape)[on_axis]
            potential[:, :count] = value
            if field:
                r0 = np.sqrt(r2+t*t)
                r1 = np.sqrt(r2+(t-length)**2)
                axial = (1/r1-1/r0)/length
                radial = (t/r0+(length-t)/r1)/(length*np.where(on_axis, 1., r2))
                electric[:, :count] = axial[..., None]*u+radial[..., None]*perp
            nodes, weights = leggauss(self.quadrature)
            for z, w in zip((nodes+1)/2, weights/2):
                y = self.starts+z*delta
                image_phi, image_e = self.images(x, y, field)
                potential[:, :count] += w*image_phi
                if field:
                    electric[:, :count] += w*image_e
        if len(self.points):
            y = self.points
            dx = x[:, None, :]-y
            radius = np.linalg.norm(dx, axis=-1)
            if np.any(radius == 0):
                raise ValueError('Evaluation point lies on an interior source point')
            potential[:, count:] = 1/radius
            if field:
                electric[:, count:] = dx/radius[..., None]**3
            image_phi, image_e = self.images(x, y, field)
            potential[:, count:] += image_phi
            if field:
                electric[:, count:] += image_e
        return (potential, electric) if field else potential


@dataclass
class VacuumSolution:
    sources: Sources
    strengths: np.ndarray
    voltage: float
    fit: dict

    def evaluate(self, points, field=False):
        points = np.asarray(points, float).reshape(-1, 3)
        phi = np.empty(len(points))
        electric = np.empty_like(points) if field else None
        for i in range(0, len(points), 128):
            kernel = self.sources.kernel(points[i:i+128], field)
            if field:
                phi[i:i+128] = kernel[0]@self.strengths
                electric[i:i+128] = np.einsum('ijk,j->ik', kernel[1], self.strengths)
            else:
                phi[i:i+128] = kernel@self.strengths
        return (phi, electric) if field else phi

    @property
    def capacitance(self):
        return 4*np.pi*epsilon_0*np.sum(self.strengths)/self.voltage


def fit_conductor(sources, boundary, weights=None, voltage=1.):
    """Unregularized least-squares Dirichlet fit; retain residual and rank."""
    if not np.isfinite(voltage) or voltage == 0:
        raise ValueError('Solve a finite nonzero bias basis; its field can then be scaled to zero')
    matrix = np.vstack([sources.kernel(boundary[i:i+128]) for i in range(0, len(boundary), 128)])
    w = np.ones(len(boundary)) if weights is None else np.sqrt(weights)
    strengths, _, rank, _ = lstsq(matrix*w[:, None], voltage*w, lapack_driver='gelsy')
    residual = matrix@strengths-voltage
    fit = dict(sources=sources.count, boundary_points=len(boundary), rank=int(rank),
               relative_rms=float(np.sqrt(np.average(residual**2, weights=w*w))/abs(voltage)),
               relative_max=float(np.max(np.abs(residual))/abs(voltage)))
    return VacuumSolution(sources, strengths, voltage, fit)


def cage_discretization(cage, panels=3, source_angles=4, axial_shift=.5, holdout=False):
    """Rounded-cylinder union boundary; holdout points are not used in fitting."""
    if panels < 1 or int(panels) != panels or source_angles < 1:
        raise ValueError('Require positive integral panel/source counts')
    starts, ends, surface, weights = [], [], [], []
    surface_angles = 11 if holdout else 12
    surface_axial = 3*panels+1 if holdout else 3*panels
    for ia, ib in cage.edges:
        a, b = cage.vertices[ia], cage.vertices[ib]
        delta = b-a
        length = np.linalg.norm(delta)
        u = delta/length
        radial = (a+b)/2
        radial /= np.linalg.norm(radial)
        other = np.cross(u, radial)
        for theta in np.arange(source_angles)*2*np.pi/source_angles:
            offset = .4*cage.wire_radius*(np.cos(theta)*radial+np.sin(theta)*other)
            for k in range(panels):
                starts.append(a+delta*k/panels+offset)
                ends.append(a+delta*(k+1)/panels+offset)
        for k in range(surface_axial):
            centre = a+delta*(k+axial_shift)/surface_axial
            angle = (np.arange(surface_angles)+(.31 if holdout else 0))*2*np.pi/surface_angles
            pts = centre+cage.wire_radius*(np.cos(angle)[:, None]*radial+np.sin(angle)[:, None]*other)
            surface.extend(pts)
            weights.extend([2*np.pi*cage.wire_radius*length/(surface_axial*surface_angles)]*len(pts))
    if holdout:
        cap_dirs = directions(181, .43)
    else:
        vertex_dirs = cage.vertices/np.linalg.norm(cage.vertices, axis=1)[:, None]
        mid = cage.vertices[cage.edges].mean(axis=1)
        cap_dirs = np.vstack((cage.cap_directions, vertex_dirs, mid/np.linalg.norm(mid, axis=1)[:, None]))
    junction_sources = []
    for i, vertex in enumerate(cage.vertices):
        pts = vertex+cage.wire_radius*cap_dirs
        surface.extend(pts)
        weights.extend([4*np.pi*cage.wire_radius**2/len(cap_dirs)]*len(pts))
        neighbours = [b if a == i else a for a, b in cage.edges if a == i or b == i]
        arms = cage.vertices[neighbours]-vertex
        arms /= np.linalg.norm(arms, axis=1)[:, None]
        junction_sources.extend([vertex, vertex+.4*cage.wire_radius*vertex/np.linalg.norm(vertex)])
        junction_sources.extend(vertex+.4*cage.wire_radius*arms)
    surface, weights = np.asarray(surface), np.asarray(weights)
    keep = cage.wire_distance(surface) >= cage.wire_radius*(1-1e-8)
    sources = Sources(np.asarray(starts), np.asarray(ends), np.asarray(junction_sources), cage.wall_radius)
    return sources, surface[keep], weights[keep]


def solve_cage(cage, panels=3, voltage=1.):
    sources, boundary, weights = cage_discretization(cage, panels)
    return fit_conductor(sources, boundary, weights, voltage)
