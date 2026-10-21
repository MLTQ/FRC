"""Limited, axis-regular reconstruction of annular volume averages."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_state import NonphysicalState
from .mhd_positivity import positive_scale, positive_blend, limiter_statistics, neighbour_minimum

_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(4)


def minmod(*values):
    same = np.logical_and.reduce([np.sign(v) == np.sign(values[0]) for v in values[1:]])
    return np.where(same, np.sign(values[0])*np.minimum.reduce([np.abs(v) for v in values]), 0)


def radial_slopes(values, coordinates):
    """MC slopes in a volume coordinate, with one-sided boundary extrapolation."""
    delta = np.diff(values, axis=0)/np.diff(coordinates, axis=0)
    out = np.empty_like(values)
    spacing = np.diff(coordinates, axis=0)
    boundary_low = delta[0]-spacing[0]/(spacing[0]+spacing[1])*(delta[1]-delta[0])
    boundary_high = delta[-1]+spacing[-1]/(spacing[-1]+spacing[-2])*(delta[-1]-delta[-2])
    out[0] = minmod(2*delta[0], boundary_low)
    out[-1] = minmod(2*delta[-1], boundary_high)
    hleft = coordinates[1:-1]-coordinates[:-2]
    hright = coordinates[2:]-coordinates[1:-1]
    centered = (hright*delta[:-1]+hleft*delta[1:])/(hleft+hright)
    out[1:-1] = minmod(2*delta[:-1], centered, 2*delta[1:])
    return out


def axial_increments(values):
    """Periodic MC slope multiplied by dz."""
    left, right = values-np.roll(values, 1, axis=1), np.roll(values, -1, axis=1)-values
    return minmod(2*left, .5*(left+right), 2*right)


def radial_pairs(low, high, *, reflect=False):
    left = np.concatenate((low[:1], high), axis=0)
    right = np.concatenate((low, high[-1:]), axis=0)
    if reflect:
        left[0, :, [1, 4]] *= -1
        right[-1, :, [1, 4]] *= -1
    return left, right


class Reconstruction:
    """Even q=a+b r²/2; odd q=r(a+b r²/2), with exact annular moments.

    Pressure slot stores total pressure during reconstruction. This retains
    constant-total-pressure straight-field equilibria without a force repair.
    No cell energy, density or pressure is changed by reconstruction.
    """
    def __init__(self, mesh, q, *, positivity=True, local_bounds=False):
        self.mesh = mesh
        # 2.5D: q has two more even, axis-regular columns, Omega = v_theta/r and g = B_theta/r.
        self.toroidal = q.shape[-1] == 8
        components = q.shape[-1]
        self.values = q.copy()
        self.values[..., 3] += (q[..., 4]**2+q[..., 5]**2)/(2*MU0)
        if self.toroidal:   # cell toroidal magnetic pressure <r^2> g^2/(2 mu0) joins the total pressure
            self.values[..., 3] += .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]*q[..., 7]**2/(2*MU0)
        self.scale = np.ones((mesh.nr, 1, components))
        self.scale[:, 0, [1, 4]] = mesh.radial_centroid[:, None]
        self.values /= self.scale
        lo, hi = mesh.rf[:-1], mesh.rf[1:]
        self.coordinates = np.broadcast_to(((lo*lo+hi*hi)/4)[:, None, None], (mesh.nr, 1, components)).copy()
        self.coordinates[:, 0, [1, 4]] = (.3*(hi**5-lo**5)/(hi**3-lo**3))[:, None]
        self.radial = radial_slopes(self.values, self.coordinates)
        if self.toroidal:
            self._wall_slopes()
        self.axial = axial_increments(self.values)
        xcenter = self.coordinates[..., 0]
        self.pressure = q[..., 3].copy()
        self.pressure_radial = radial_slopes(self.pressure, xcenter)
        self.pressure_axial = axial_increments(self.pressure)
        density_scale = np.ones_like(self.pressure)
        thermal_scale = np.ones_like(self.pressure)
        # Optional local bounds: density/thermal-pressure traces stay >= the 3x3 neighbour minimum,
        # so trace wave speeds cannot exceed the cell speeds that set the time step by a large factor.
        density_bound = neighbour_minimum(q[..., 0]) if positivity and local_bounds else None
        pressure_bound = neighbour_minimum(q[..., 3]) if positivity and local_bounds else None
        if positivity:
            density_scale = positive_scale(q[..., 0], self._linear_minimum(q[..., 0], self.radial[..., 0], self.axial[..., 0]),
                                           density_bound)
            self.radial[..., 0] *= density_scale
            self.axial[..., 0] *= density_scale
            thermal_scale = positive_scale(self.pressure, self._linear_minimum(self.pressure, self.pressure_radial, self.pressure_axial),
                                           pressure_bound)
            self.pressure_radial *= thermal_scale
            self.pressure_axial *= thermal_scale
        self.pressure_weight = np.ones_like(self.pressure)
        # Each trace carries the radius at which its toroidal magnetic pressure (R g)^2/(2 mu0) is evaluated:
        # the face radius radially, sqrt(<r^2>) on axial faces (as in toroidal_flux).
        self.axial_radius = np.sqrt(.5*(lo*lo+hi*hi))[:, None]
        radial = {(upper, offset): (self._raw_radial(mesh.rf[1:] if upper else mesh.rf[:-1], offset),
                                     self._auxiliary_pressure(mesh.rf[1:] if upper else mesh.rf[:-1], offset),
                                     (mesh.rf[1:] if upper else mesh.rf[:-1])[:, None])
                  for upper in (False, True) for offset in (-.5, 0., .5)}
        axial = {offset: ((self.values+offset*self.axial)*self.scale, self.pressure+offset*self.pressure_axial,
                          self.axial_radius)
                 for offset in (-.5, .5)}
        if positivity:
            for raw, auxiliary, radius in (*radial.values(), *axial.values()):
                proposed = raw[..., 3]-(raw[..., 4]**2+raw[..., 5]**2)/(2*MU0)
                if self.toroidal:
                    proposed = proposed-(radius*raw[..., 7])**2/(2*MU0)
                self.pressure_weight = np.minimum(self.pressure_weight, positive_blend(auxiliary, proposed, self.pressure,
                                                                                         pressure_bound))
        self.statistics = limiter_statistics(density_scale, thermal_scale, self.pressure_weight)
        self._radial_states = {key: self.primitive(*value) for key, value in radial.items()}
        self._axial_states = {key: self.primitive(*value) for key, value in axial.items()}

    def _wall_slopes(self):
        """Omega at a gas wall: MC slope against the mirrored (even: slip) ghost, which is identically zero,
        so the wall cell's Omega is piecewise constant (first order there).

        The default one-sided extrapolation allows twice the interior difference. Nothing restores a
        wall-cell extremum of the angular velocity, and with the large fast speed of the tenuous gap the
        Rusanov dissipation at the inner face then pumps angular momentum into the wall cell (measured in
        the first Hall formation run: v_theta 2e6 m/s in the 1%-of-fill gap, while the neighbouring cells
        turned the other way). The ghost-bounded slope keeps the inner-face trace between the two cell
        means. g keeps the one-sided slope: bounding it too made the torsional Alfven wave gain 3% of its
        energy per period at 16 cells. The axis (regular, x = r^2/2 >= 0) is unchanged."""
        self.radial[-1, :, 6] = 0.
        if self.mesh.annular:
            self.radial[0, :, 6] = 0.

    def _linear_minimum(self, center, radial, axial):
        x = self.coordinates[..., 0]
        lo = radial*(.5*self.mesh.rf[:-1, None]**2-x)
        hi = radial*(.5*self.mesh.rf[1:, None]**2-x)
        return center+np.minimum(lo, hi)-.5*np.abs(axial)

    def _auxiliary_pressure(self, radius, offset=0):
        return self.pressure+self.pressure_radial*(.5*radius[:, None]**2-self.coordinates[..., 0])+offset*self.pressure_axial

    def primitive(self, values, auxiliary, radius=None):
        out = values.copy()
        out[..., 3] -= (out[..., 4]**2+out[..., 5]**2)/(2*MU0)
        if self.toroidal:
            out[..., 3] -= (radius*out[..., 7])**2/(2*MU0)
        limited = auxiliary+self.pressure_weight*(out[..., 3]-auxiliary)
        out[..., 3] = np.where(self.pressure_weight < 1, limited, out[..., 3])
        if not np.all(np.isfinite(out)) or np.any(out[..., 0] <= 0) or np.any(out[..., 3] <= 0):
            raise NonphysicalState(f"Reconstructed face/corner state invalid: minimum density {out[..., 0].min():.6g} kg/m^3, thermal pressure {out[..., 3].min():.6g} Pa; no floor applied")
        return out

    def _raw_radial(self, radius, axial_offset=0):
        values = self.values+self.radial*(.5*radius[:, None, None]**2-self.coordinates)+axial_offset*self.axial
        factor = np.ones_like(self.scale)
        factor[:, 0, [1, 4]] = radius[:, None]
        return values*factor

    def radial_values(self, upper, axial_offset=0):
        return self._radial_states[(upper, axial_offset)]

    def faces(self):
        left, right = radial_pairs(self.radial_values(False), self.radial_values(True), reflect=True)
        low, high = self._axial_states[-.5], self._axial_states[.5]
        down = np.concatenate((high[:, -1:], high), axis=1)
        up = np.concatenate((low, low[:, :1]), axis=1)
        return left, right, down, up

    def corners(self):
        left_down, right_down = radial_pairs(self.radial_values(False, -.5), self.radial_values(True, -.5), reflect=True)
        left_up, right_up = radial_pairs(self.radial_values(False, .5), self.radial_values(True, .5), reflect=True)
        sw = np.concatenate((left_up[:, -1:], left_up), axis=1)
        se = np.concatenate((right_up[:, -1:], right_up), axis=1)
        nw = np.concatenate((left_down, left_down[:, :1]), axis=1)
        ne = np.concatenate((right_down, right_down[:, :1]), axis=1)
        return sw, se, nw, ne

    def radial_source(self):
        # Integrate pt dr / integral(r dr), NOT pt/r at an assumed point.
        original = (self.values[..., 3]-self.radial[..., 3]*self.mesh.dr**2/12)/self.mesh.r[:, None]
        if np.all(self.pressure_weight == 1):
            return original+self._toroidal_source() if self.toroidal else original
        auxiliary = (self.pressure-self.pressure_radial*self.mesh.dr**2/12)/self.mesh.r[:, None]
        # Br is cubic and Bz quadratic in r: four-point Gauss integrates B² dr
        # exactly. Use the same blended pressure as the face/corner states.
        magnetic = np.zeros_like(self.pressure)
        for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS):
            radius = self.mesh.r+.5*self.mesh.dr*node
            raw = self._raw_radial(radius)
            magnetic += weight*(raw[..., 4]**2+raw[..., 5]**2)/(2*MU0)
            if self.toroidal:
                magnetic += weight*(radius[:, None]*raw[..., 7])**2/(2*MU0)
        auxiliary += .5*self.mesh.dr*magnetic/self.mesh.radial_measure[:, None]
        blended = self.pressure_weight*original+(1-self.pressure_weight)*auxiliary
        return blended+self._toroidal_source() if self.toroidal else blended

    def _toroidal_source(self):
        """Centrifugal and hoop terms, integral of (rho v_theta^2 - B_theta^2/mu0) dr / integral(r dr), by
        four-point Gauss on the reconstructed rho, Omega and g (v_theta = r Omega, B_theta = r g)."""
        total = np.zeros_like(self.pressure)
        for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS):
            radius = (self.mesh.r+.5*self.mesh.dr*node)[:, None]
            raw = self._raw_radial(radius[:, 0])
            total += weight*radius**2*(raw[..., 0]*raw[..., 6]**2-raw[..., 7]**2/MU0)
        return .5*self.mesh.dr*total/self.mesh.radial_measure[:, None]


def staggered_corner_fields(mesh, br_faces, bz_faces):
    """Reconstruct each area-averaged field only transverse to its face normal."""
    x = ((mesh.rf[:-1]**2+mesh.rf[1:]**2)/4)[:, None]
    slope = radial_slopes(bz_faces, x)
    low = bz_faces+slope*(.5*mesh.rf[:-1, None]**2-x)
    high = bz_faces+slope*(.5*mesh.rf[1:, None]**2-x)
    bz_left, bz_right = radial_pairs(low, high)
    increment = axial_increments(br_faces)
    br_low, br_high = br_faces-.5*increment, br_faces+.5*increment
    br_down = np.concatenate((br_high[:, -1:], br_high), axis=1)
    br_up = np.concatenate((br_low, br_low[:, :1]), axis=1)
    return bz_left, bz_right, br_down, br_up
