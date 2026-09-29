"""Regular axisymmetric radial FE / axial Fourier linear Hall cylinder."""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.constants import mu_0 as MU0
from .hall_modes import LinearHallModes
from .hall_parameters import HallParameters


@dataclass(frozen=True)
class HallCylinderMesh:
    radius: float
    length: float
    nr: int
    nz: int
    quadrature: int = 8

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in (self.radius, self.length)):
            raise ValueError('Positive finite cylinder dimensions required')
        if any(not isinstance(n, int) or isinstance(n, bool) or n < 4 for n in (self.nr, self.nz)):
            raise ValueError('At least four radial cells and five axial samples required')
        if self.nz % 2 == 0:
            raise ValueError('Odd axial count required to avoid an ambiguous Nyquist mode')
        if not isinstance(self.quadrature, int) or isinstance(self.quadrature, bool) or self.quadrature < 2:
            raise ValueError('At least two radial quadrature points required')

    @property
    def r(self):
        return np.linspace(0., self.radius, self.nr+1)

    @property
    def z(self):
        return np.arange(self.nz)*self.length/self.nz


class LinearHallCylinder:
    def __init__(self, mesh: HallCylinderMesh, parameters: HallParameters):
        self.mesh, self.parameters = mesh, parameters
        r = mesh.r
        g, weight = leggauss(mesh.quadrature)
        mass = np.zeros((mesh.nr+1, mesh.nr+1))
        stiffness = np.zeros_like(mass)
        for cell, (left, right) in enumerate(zip(r[:-1], r[1:])):
            points = left+(g+1)*(right-left)/2
            weights = weight*(right-left)/2*points
            shape = np.array([(1-g)/2, (1+g)/2]).T
            derivative = np.broadcast_to(np.array([-1., 1.])/(right-left), shape.shape)
            cylindrical = derivative+shape/points[:, None]
            index = np.ix_([cell, cell+1], [cell, cell+1])
            mass[index] += np.einsum('q,qi,qj->ij', weights, shape, shape)
            stiffness[index] += np.einsum('q,qi,qj->ij', weights, cylindrical, cylindrical)
        self.mass = mass[1:-1, 1:-1]
        self.stiffness = stiffness[1:-1, 1:-1]
        self.radial_eigenvalues, self.basis = eigh(self.stiffness, self.mass)
        self.kz = 2*np.pi*np.fft.fftfreq(mesh.nz, d=mesh.length/mesh.nz)
        self.modes = LinearHallModes(self.kz, self.radial_eigenvalues, parameters)
        self.root = np.sqrt(self.modes.k_squared)
        self.volume_factor = 2*np.pi*mesh.length
        self.velocity_scale = np.sqrt(MU0*parameters.density)

    def from_nodal(self, a, toroidal, stream, azimuthal_velocity):
        """Interpolate resolved real nodal potentials/fields; endpoints must be zero."""
        fields = []
        for field in (a, toroidal, stream, azimuthal_velocity):
            field = np.asarray(field)
            if field.shape != (self.mesh.nr+1, self.mesh.nz) or not np.isrealobj(field) or not np.isfinite(field).all():
                raise ValueError('Finite real nodal fields with shape (nr+1,nz) required')
            if np.any(field[[0, -1]] != 0):
                raise ValueError('Nodal potential/toroidal traces must vanish at axis and outer boundary')
            fields.append(self.basis.T@self.mass@np.fft.fft(field[1:-1], axis=1, norm='forward'))
        fields[0] *= self.root
        fields[2] *= self.root*self.velocity_scale
        fields[3] *= self.velocity_scale
        return np.stack(fields, axis=-1)

    def _validate(self, state):
        state = self.modes.validate(state)
        partner = (-np.arange(self.mesh.nz)) % self.mesh.nz
        error = np.max(abs(state-state[:, partner].conj()))
        scale = np.max(abs(state))
        if error > 256*np.finfo(float).eps*scale:
            raise ValueError('Cylinder state must have the conjugate symmetry of a real field')
        return state

    def step(self, state, dt):
        self._validate(state)
        return self.modes.step(state, dt, self.volume_factor)

    def energies(self, state):
        self._validate(state)
        return self.modes.energies(state, self.volume_factor)

    def sample(self, state, r, z, quantity='magnetic'):
        """Evaluate real physical fields at broadcast-compatible points, SI units."""
        state = self._validate(state)
        r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
        if not np.isfinite(r).all() or not np.isfinite(z).all() or np.any((r < 0) | (r > self.mesh.radius)):
            raise ValueError('Finite sample points inside cylinder radius required')
        output_shape = r.shape+(3,)
        r, z = r.ravel(), z.ravel()
        cell = np.minimum((r/(self.mesh.radius/self.mesh.nr)).astype(int), self.mesh.nr-1)
        fraction = (r-self.mesh.r[cell])/(self.mesh.radius/self.mesh.nr)
        extended = np.pad(self.basis, ((1, 1), (0, 0)))
        values = (1-fraction[:, None])*extended[cell]+fraction[:, None]*extended[cell+1]
        derivatives = (extended[cell+1]-extended[cell])/(self.mesh.radius/self.mesh.nr)
        if quantity == 'magnetic':
            potential, toroidal = state[..., 0]/self.root, state[..., 1]
        elif quantity == 'velocity':
            potential = state[..., 2]/self.root/self.velocity_scale
            toroidal = state[..., 3]/self.velocity_scale
        elif quantity == 'current':
            potential, toroidal = state[..., 1]/MU0, self.root*state[..., 0]/MU0
        else:
            raise ValueError('quantity must be magnetic, velocity, or current')
        phase = np.exp(1j*z[:, None]*self.kz)
        def evaluate(radial, coefficients):
            return np.einsum('pm,mk,pk->p', radial, coefficients, phase).real
        axial_basis = derivatives+np.divide(values, r[:, None], out=derivatives.copy(), where=r[:, None] != 0)
        radial = evaluate(values, -1j*self.kz*potential)
        theta = evaluate(values, toroidal)
        axial = evaluate(axial_basis, potential)
        return np.stack((radial, theta, axial), axis=-1).reshape(output_shape)

    def quadrature_points(self, order=8, axial_count=None):
        """Independent caller-selectable physical volume quadrature, axis excluded."""
        if not isinstance(order, int) or isinstance(order, bool) or order < 2:
            raise ValueError('Positive integration order at least two required')
        axial_count = self.mesh.nz if axial_count is None else axial_count
        if not isinstance(axial_count, int) or axial_count < self.mesh.nz:
            raise ValueError('Axial quadrature must resolve at least the stored spectrum')
        x, w = leggauss(order)
        dr = self.mesh.radius/self.mesh.nr
        r = (self.mesh.r[:-1, None]+(x+1)*dr/2).ravel()
        weights = (np.broadcast_to(w, (self.mesh.nr, order))*dr/2).ravel()*r
        z = np.arange(axial_count)*self.mesh.length/axial_count
        return r[:, None], z[None, :], 2*np.pi*weights[:, None]*self.mesh.length/axial_count*np.ones((1, axial_count))

    def angular_momentum(self, state):
        r, z, weights = self.quadrature_points(order=2)
        return float(np.sum(weights*self.parameters.density*r*self.sample(state, r, z, 'velocity')[..., 1]))
