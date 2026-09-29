"""Conservative cylindrical charge fluxes dual to toroidal Faraday work."""
import numpy as np
from scipy.constants import mu_0 as MU0


class DualCurrent:
    """Nodal dual volumes, periodic z; a physical toroidal wall trace is explicit."""
    def __init__(self, mesh):
        self.mesh = mesh
        self.x = .5*mesh.rf**2
        self.dx = np.diff(self.x)
        self.dual_x_faces = np.r_[self.x[0], .5*(self.x[:-1]+self.x[1:]), self.x[-1]]
        self.volumes = np.broadcast_to(2*np.pi*np.diff(self.dual_x_faces)[:, None]*mesh.dz, (mesh.nr+1, mesh.nz)).copy()
        self.r, self.z = np.meshgrid(mesh.rf, mesh.zf[:-1], indexing='ij')
        self.magnetic_mass = mesh.volumes*(self.x[:-1]+self.x[1:])[:, None]/MU0

    def scalar(self, values, name='dual scalar'):
        result = np.asarray(values, float)
        if result.shape != self.volumes.shape or not np.isfinite(result).all():
            raise ValueError(f'Finite {name} on all periodic dual nodes required')
        return result

    def current_flux(self, g, wall_btheta):
        """Return integrated current [A] through radial/axial dual faces."""
        mesh = self.mesh
        g = np.asarray(g, float)
        if g.shape != (mesh.nr, mesh.nz) or not np.isfinite(g).all():
            raise ValueError('Finite primal cell g required')
        wall = np.broadcast_to(np.asarray(wall_btheta, float), (mesh.nz,))
        if not np.isfinite(wall).all():
            raise ValueError('Finite physical wall Btheta trace required')
        stream = 2*np.pi*(self.x[:-1]+self.x[1:])[:, None]*g/MU0
        boundary = 2*np.pi*mesh.radius*wall/MU0
        radial = np.zeros((mesh.nr+2, mesh.nz))
        radial[1:-1] = np.roll(stream, 1, axis=1)-stream
        radial[-1] = np.roll(boundary, 1)-boundary
        axial = np.diff(np.concatenate((np.zeros((1, mesh.nz)), stream, boundary[None]), axis=0), axis=0)
        return radial, axial

    def divergence(self, radial, axial):
        """Outward integrated face divergence, divided by the dual volume."""
        return (np.diff(radial, axis=0)+axial-np.roll(axial, 1, axis=1))/self.volumes

    def gradient(self, scalar):
        """Er/r and Ez coefficients of a bilinear scalar in x=r²/2,z."""
        scalar = self.scalar(scalar)
        return np.diff(scalar, axis=0)/self.dx[:, None], (np.roll(scalar, -1, axis=1)-scalar)/self.mesh.dz

    def faraday(self, radial_coefficient, axial_electric):
        """Physical Btheta/r rate; both radial endpoint Ez values are retained."""
        return np.diff(axial_electric, axis=0)/self.dx[:, None]-(np.roll(radial_coefficient, -1, axis=1)-radial_coefficient)/self.mesh.dz

    def work(self, g, wall_btheta, radial_coefficient, axial_electric):
        radial, axial = self.current_flux(g, wall_btheta)
        volume = float(np.sum(radial[1:-1]*self.dx[:, None]*radial_coefficient)+np.sum(axial*self.mesh.dz*axial_electric))
        boundary = float(np.sum(2*np.pi*self.mesh.radius*np.asarray(wall_btheta)*self.mesh.dz*axial_electric[-1]/MU0))
        magnetic = float(np.sum(self.magnetic_mass*g*self.faraday(radial_coefficient, axial_electric)))
        return dict(current_work_w=volume, boundary_magnetic_input_w=boundary, magnetic_rate_w=magnetic)
