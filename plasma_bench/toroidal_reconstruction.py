"""Axis-regular reconstruction of rotation, toroidal field and total pressure."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_reconstruction import radial_slopes, axial_increments, radial_pairs
from .mhd_positivity import positive_scale, positive_blend, limiter_statistics
from .mhd_state import NonphysicalState
from .toroidal_state import radius_squared

NODES, WEIGHTS = np.polynomial.legendre.leggauss(5)


class ToroidalReconstruction:
    def __init__(self, mesh, q, *, order=2, positivity=True):
        self.mesh = mesh
        self.values = q.copy()
        self.values[..., 3] += (q[..., 4]**2+q[..., 5]**2+radius_squared(mesh)*q[..., 7]**2)/(2*MU0)
        self.scale = np.ones((mesh.nr, 1, 8))
        lo, hi = mesh.rf[:-1], mesh.rf[1:]
        self.coordinates = np.broadcast_to(((lo**2+hi**2)/4)[:, None, None], (mesh.nr, 1, 8)).copy()
        if order == 2:
            self.scale[:, 0, [1, 4]] = mesh.radial_centroid[:, None]
            self.coordinates[:, 0, [1, 4]] = (.3*(hi**5-lo**5)/(hi**3-lo**3))[:, None]
        self.regular_poloidal = order == 2
        self.values /= self.scale
        self.radial = radial_slopes(self.values, self.coordinates) if order == 2 else np.zeros_like(self.values)
        self.axial = axial_increments(self.values) if order == 2 else np.zeros_like(self.values)
        self.pressure = q[..., 3].copy()
        self.pressure_radial = radial_slopes(self.pressure, self.coordinates[..., 0]) if order == 2 else np.zeros_like(self.pressure)
        self.pressure_axial = axial_increments(self.pressure) if order == 2 else np.zeros_like(self.pressure)
        density_scale = np.ones_like(self.pressure)
        thermal_scale = np.ones_like(self.pressure)
        if positivity:
            density_scale = positive_scale(q[..., 0], self.minimum(q[..., 0], self.radial[..., 0], self.axial[..., 0]))
            self.radial[..., 0] *= density_scale
            self.axial[..., 0] *= density_scale
            thermal_scale = positive_scale(self.pressure, self.minimum(self.pressure, self.pressure_radial, self.pressure_axial))
            self.pressure_radial *= thermal_scale
            self.pressure_axial *= thermal_scale
        self.pressure_weight = np.ones_like(self.pressure)
        radial = {(upper, offset): (self.raw(mesh.rf[1:] if upper else mesh.rf[:-1], offset),
                                    self.auxiliary(mesh.rf[1:] if upper else mesh.rf[:-1], offset),
                                    (mesh.rf[1:] if upper else mesh.rf[:-1])[:, None])
                  for upper in (False, True) for offset in (-.5, 0., .5)}
        axial = {offset: ((self.values+offset*self.axial)*self.scale, self.pressure+offset*self.pressure_axial,
                          np.sqrt(radius_squared(mesh))) for offset in (-.5, .5)}
        gauss = [(self.raw(mesh.r+.5*mesh.dr*x), self.auxiliary(mesh.r+.5*mesh.dr*x),
                  (mesh.r+.5*mesh.dr*x)[:, None]) for x in NODES]
        if positivity:
            for values, auxiliary, r in (*radial.values(), *axial.values(), *gauss):
                p = values[..., 3]-self.magnetic_pressure(values, r)
                self.pressure_weight = np.minimum(self.pressure_weight, positive_blend(auxiliary, p, self.pressure))
        self.statistics = limiter_statistics(density_scale, thermal_scale, self.pressure_weight)
        self.radial_states = {key: self.primitive(*args) for key, args in radial.items()}
        self.axial_states = {key: self.primitive(*args) for key, args in axial.items()}
        self.gauss_states = [(self.primitive(*args), args[2]) for args in gauss]

    @staticmethod
    def magnetic_pressure(values, r):
        return (values[..., 4]**2+values[..., 5]**2+(r*values[..., 7])**2)/(2*MU0)

    def minimum(self, center, radial, axial):
        x = self.coordinates[..., 0]
        return center+np.minimum(radial*(.5*self.mesh.rf[:-1, None]**2-x),
                                 radial*(.5*self.mesh.rf[1:, None]**2-x))-.5*abs(axial)

    def auxiliary(self, radius, offset=0.):
        return self.pressure+self.pressure_radial*(.5*radius[:, None]**2-self.coordinates[..., 0])+offset*self.pressure_axial

    def raw(self, radius, offset=0.):
        values = self.values+self.radial*(.5*radius[:, None, None]**2-self.coordinates)+offset*self.axial
        if self.regular_poloidal:
            factor = np.ones_like(self.scale)
            factor[:, 0, [1, 4]] = radius[:, None]
            values = values*factor
        return values

    def primitive(self, values, auxiliary, r):
        out = values.copy()
        out[..., 3] -= self.magnetic_pressure(out, r)
        limited = auxiliary+self.pressure_weight*(out[..., 3]-auxiliary)
        out[..., 3] = np.where(self.pressure_weight < 1, limited, out[..., 3])
        if not np.isfinite(out).all() or np.any(out[..., 0] <= 0) or np.any(out[..., 3] <= 0):
            raise NonphysicalState('Invalid toroidal face/corner/source trace; no floor applied')
        return out

    def faces(self):
        left, right = radial_pairs(self.radial_states[False, 0.], self.radial_states[True, 0.], reflect=True)
        low, high = self.axial_states[-.5], self.axial_states[.5]
        down = np.concatenate((high[:, -1:], high), axis=1)
        up = np.concatenate((low, low[:, :1]), axis=1)
        return left, right, down, up

    def corners(self):
        ld, rd = radial_pairs(self.radial_states[False, -.5], self.radial_states[True, -.5], reflect=True)
        lu, ru = radial_pairs(self.radial_states[False, .5], self.radial_states[True, .5], reflect=True)
        return (np.concatenate((lu[:, -1:], lu), axis=1), np.concatenate((ru[:, -1:], ru), axis=1),
                np.concatenate((ld, ld[:, :1]), axis=1), np.concatenate((rd, rd[:, :1]), axis=1))

    def radial_source(self):
        source = np.zeros_like(self.pressure)
        for weight, (q, r) in zip(WEIGHTS, self.gauss_states):
            source += weight*(q[..., 3]+self.magnetic_pressure(q, r)+q[..., 0]*(r*q[..., 6])**2-(r*q[..., 7])**2/MU0)
        return .5*self.mesh.dr*source/self.mesh.radial_measure[:, None]
