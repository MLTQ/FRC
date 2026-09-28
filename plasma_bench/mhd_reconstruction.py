"""Limited, axis-regular reconstruction of annular volume averages."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_state import NonphysicalState


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
    def __init__(self, mesh, q):
        self.mesh = mesh
        self.values = q.copy()
        self.values[..., 3] += (q[..., 4]**2+q[..., 5]**2)/(2*MU0)
        self.scale = np.ones((mesh.nr, 1, 6))
        self.scale[:, 0, [1, 4]] = mesh.radial_centroid[:, None]
        self.values /= self.scale
        lo, hi = mesh.rf[:-1], mesh.rf[1:]
        self.coordinates = np.broadcast_to(((lo*lo+hi*hi)/4)[:, None, None], (mesh.nr, 1, 6)).copy()
        self.coordinates[:, 0, [1, 4]] = (.3*(hi**5-lo**5)/(hi**3-lo**3))[:, None]
        self.radial = radial_slopes(self.values, self.coordinates)
        self.axial = axial_increments(self.values)

    @staticmethod
    def primitive(values):
        out = values.copy()
        out[..., 3] -= (out[..., 4]**2+out[..., 5]**2)/(2*MU0)
        if not np.all(np.isfinite(out)) or np.any(out[..., 0] <= 0) or np.any(out[..., 3] <= 0):
            raise NonphysicalState(f"Reconstructed face/corner state invalid: minimum density {out[..., 0].min():.6g} kg/m^3, thermal pressure {out[..., 3].min():.6g} Pa; no floor applied")
        return out

    def radial_values(self, upper, axial_offset=0):
        radius = self.mesh.rf[1:] if upper else self.mesh.rf[:-1]
        values = self.values+self.radial*(.5*radius[:, None, None]**2-self.coordinates)+axial_offset*self.axial
        factor = np.ones_like(self.scale)
        factor[:, 0, [1, 4]] = radius[:, None]
        return values*factor

    def faces(self):
        left, right = radial_pairs(self.radial_values(False), self.radial_values(True), reflect=True)
        low = (self.values-.5*self.axial)*self.scale
        high = (self.values+.5*self.axial)*self.scale
        down = np.concatenate((high[:, -1:], high), axis=1)
        up = np.concatenate((low, low[:, :1]), axis=1)
        return tuple(self.primitive(v) for v in (left, right, down, up))

    def corners(self):
        left_down, right_down = radial_pairs(self.radial_values(False, -.5), self.radial_values(True, -.5), reflect=True)
        left_up, right_up = radial_pairs(self.radial_values(False, .5), self.radial_values(True, .5), reflect=True)
        sw = np.concatenate((left_up[:, -1:], left_up), axis=1)
        se = np.concatenate((right_up[:, -1:], right_up), axis=1)
        nw = np.concatenate((left_down, left_down[:, :1]), axis=1)
        ne = np.concatenate((right_down, right_down[:, :1]), axis=1)
        return tuple(self.primitive(v) for v in (sw, se, nw, ne))

    def radial_source(self):
        # Integrate pt dr / integral(r dr), NOT pt/r at an assumed point.
        return (self.values[..., 3]-self.radial[..., 3]*self.mesh.dr**2/12)/self.mesh.r[:, None]


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
