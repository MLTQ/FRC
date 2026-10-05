"""Small axial-Fourier symbols of the physical uniform plasma preconditioner."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.linalg import block_diag
from scipy.sparse import diags, eye, kron
from .compatible_geometry import gram
from .regular_fluid_linear import radial_matrices, periodic_derivative


class ModalWallBlocks:
    def __init__(self, coordinates):
        self.coordinates = c = coordinates; d = c.fluid; geo = c.geometry; mesh = geo.mesh
        self.radial_mass, self.radial_gradient = radial_matrices(d)
        background = c.background(); fields = d.fields(background)
        entropy = fields['entropy'][0, 0, 0]; h = fields['specific_internal'][0, 0, 0]; v = c.axial_speed
        pressure = np.zeros(6); pressure[0] = np.sum((d.thermal.gammas-1)*h*(d.thermal.gammas-entropy))
        pressure[4:] = (d.thermal.gammas-1)*h
        self.kr = np.zeros((6, 6)); self.kr[[0, 3, 4, 5], 1] = [1., v, *entropy]
        self.kp = np.zeros((6, 6)); self.kp[1] = pressure/d.mass
        self.kz = np.zeros((6, 6)); self.kz[0, 3] = 1.; self.kz[1, 1] = self.kz[2, 2] = v
        self.kz[3, 3] = 2*v; self.kz[3, 0] = -v*v; self.kz[3] += pressure/d.mass
        for k in range(2):
            self.kz[4+k, 0] = -entropy[k]*v; self.kz[4+k, 3] = entropy[k]; self.kz[4+k, 4+k] = v
        z = periodic_derivative(d.nz, d.dz, d.coefficients)
        self.z_symbol = np.fft.fft(z[:, 0].toarray().ravel())
        self.node_symbol = np.fft.fft(gram(geo.Nz, geo.wz)[:, 0].toarray().ravel())
        self.edge_symbol = np.fft.fft(gram(geo.Ez, geo.wz)[:, 0].toarray().ravel())
        self.derivative_symbol = np.fft.fft(gram(geo.dNz, geo.wz)[:, 0].toarray().ravel())
        self.radial_e_mass = gram(geo.radial_e, geo.wr).toarray()
        self.theta_e_mass = gram(geo.theta_e, geo.wr).toarray()
        self.axial_e_mass = gram(geo.Nx, geo.wr).toarray()
        self.psi_radial = gram(geo.radial_b[:, 1:-1], geo.wr).toarray()
        self.psi_axial = gram(geo.dNx[:, 1:-1], geo.wr).toarray()
        self.cross_e = (geo.radial_e.T@geo.theta_e.multiply(geo.wr[:, None])).toarray()
        radial_velocity = c.system.mapping.radial_velocity
        self.force_theta = (radial_velocity.T@geo.theta_e.multiply((geo.wr*geo.radial_points)[:, None])).toarray()
        self.force_radial = (radial_velocity.T@geo.radial_e.multiply((geo.wr*geo.radial_points)[:, None])).toarray()
        self.force_denominator = (d.mass*d.volume[..., 0]*d.r[..., 0]**2).ravel()
        kinetic = np.zeros((d.elements*(d.degree+1), 6))
        kinetic[:, 1:3] = (d.volume[..., 0]*d.mass*d.r[..., 0]**2/c.density).reshape(-1, 1)
        self.kinetic = kinetic.ravel()[c.fluid_rows]*c.fluid_row_scale**2
        self.lift = block_diag(np.eye(mesh.nr-1)*c.magnetic_scalars[0], c.radial_basis*c.magnetic_scalars[1])
        self.reduced_mass = c.system.ion_mass*ME/(c.system.ion_mass+ME)
        self.delta = (c.system.ion_mass-ME)/(c.system.ion_mass+ME)

    def fluid(self, mode):
        c = self.coordinates
        a = (kron(self.radial_mass, self.kr)+kron(self.radial_gradient, self.kp)
             +kron(eye(self.radial_mass.shape[0]), self.z_symbol[mode]*self.kz)).tocsr()
        a = a[c.fluid_rows][:, c.fluid_rows]
        return (diags(1/c.fluid_row_scale)@a@diags(c.fluid_row_scale)).tocsc()

    def coupled(self, mode):
        c = self.coordinates; geo = c.geometry; mesh = geo.mesh; nr = mesh.nr
        zn, ze, zd = self.node_symbol[mode], self.edge_symbol[mode], self.derivative_symbol[mode]
        mass = block_diag(self.radial_e_mass*zn, self.theta_e_mass*zn, self.axial_e_mass*ze)
        hessian = block_diag((self.psi_radial*zd+self.psi_axial*zn)/MU0, self.radial_e_mass*ze/MU0)
        curl = np.zeros((2*nr-1, 3*nr), complex)
        curl[:nr-1, nr:2*nr-1] = -np.diag(mesh.rf[1:-1])
        curl[nr-1:, :nr] = np.eye(nr)*(1-np.exp(2j*np.pi*mode/c.nz))/mesh.dz
        curl[nr-1:, 2*nr-1:] = np.diff(np.eye(nr+1), axis=0)/mesh.radial_measure[:, None]
        current = -np.linalg.solve(mass, curl.conj().T@hessian@self.lift)
        metric = self.lift.conj().T@hessian@self.lift+self.reduced_mass/(E**2*c.density)*current.conj().T@mass@current
        hall = np.zeros_like(mass); hall[:nr, nr:2*nr-1] = self.cross_e*zn
        hall[nr:2*nr-1, :nr] = -self.cross_e.T*zn
        magnetic = np.linalg.solve(metric, -self.delta*c.bias/(E*c.density)*current.conj().T@hall@current)
        magnetic += np.eye(magnetic.shape[0])*c.axial_speed*self.z_symbol[mode]
        phase = np.exp(2j*np.pi*mode*np.arange(c.nz)/c.nz)
        fluid_shape = c.system.mapping.axial@phase; field_shape = geo.Nz@phase
        cross_z = np.sum(geo.wz*fluid_shape.conj()*field_shape)/c.nz
        force = np.zeros((len(self.force_denominator), 6, 3*nr), complex)
        active = self.force_denominator > 0
        force[active, 1, nr:2*nr-1] = c.bias*cross_z*self.force_theta[active]/self.force_denominator[active, None]
        force[active, 2, :nr] = -c.bias*cross_z*self.force_radial[active]/self.force_denominator[active, None]
        force = force.reshape(-1, 3*nr)[c.fluid_rows]/c.fluid_row_scale[:, None]
        to_fluid = force@current
        to_magnetic = np.linalg.solve(metric, -to_fluid.conj().T*self.kinetic[None, :])
        return to_fluid, to_magnetic, magnetic
