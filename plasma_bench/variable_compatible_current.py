"""Physical variable-density inertia and local current/magnetic energy rates."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.sparse import block_diag
from scipy.sparse.linalg import LinearOperator
from .compatible_geometry import CompatibleGeometry
from .compatible_density import CompatibleDensity


def weak_vector(geometry, values):
    """Physical tensor quadrature adjoint P.T W, without a mass solve."""
    values=np.asarray(values,float)
    if values.shape!=geometry.r.shape+(3,) or not np.isfinite(values).all():
        raise ValueError('Finite physical vector on common quadrature required')
    result=[]
    for k,(a,b) in enumerate(geometry.pairs):
        weighted=geometry._unlayout(geometry.weight*values[...,k])
        result.append(np.asarray((a.T@weighted)@b).ravel())
    return np.concatenate(result)


class VariableCompatibleCurrent:
    """Instantaneous component with supplied conservative density, no density step."""
    def __init__(self, mesh, number_density, ion_mass, *, degree=5,
                 quadrature_order=10, geometry=None):
        if not np.isfinite(ion_mass) or ion_mass<=0:
            raise ValueError('Positive finite physical ion mass required')
        self.mesh=mesh
        self.geometry=geometry if geometry is not None else CompatibleGeometry(
            mesh,degree=degree,quadrature_order=quadrature_order)
        if self.geometry.mesh!=mesh:raise ValueError('Matching density and field mesh required')
        self.density_map=CompatibleDensity(self.geometry,number_density)
        self.density=self.density_map.samples
        self.ion_mass=ion_mass;self.electron_mass=ME
        self.alpha=ME/(ion_mass+ME);self.beta=ion_mass/(ion_mass+ME)
        self.delta=self.beta-self.alpha;self.reduced_mass=ion_mass*self.alpha
        with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
            self.inertia_coefficient=self.reduced_mass/(self.density*E**2)
        if not np.isfinite(self.inertia_coefficient).all() or np.any(self.inertia_coefficient<=0):
            raise ValueError('Physical current inertia exceeds positive finite numerical range')
        self.inertia_coefficient.setflags(write=False)
        self.mass=block_diag(self.geometry.mass,format='csc');self.hessian=self.geometry.hessian
        self.inertia=LinearOperator(self.mass.shape,matvec=self._inertia,
                                    rmatvec=self._inertia,dtype=float)

    def _inertia(self, current):
        j=self.geometry.sample_electric(np.asarray(current).reshape(-1))
        return weak_vector(self.geometry,self.inertia_coefficient[...,None]*j)

    def current(self, magnetic):
        magnetic=np.asarray(magnetic,float)
        if magnetic.shape!=(self.hessian.shape[0],) or not np.isfinite(magnetic).all():
            raise ValueError('Finite complete magnetic increment required')
        weak=-self.geometry.faraday.T@(self.hessian@magnetic)
        return np.concatenate([solve(v) for solve,v in zip(
            self.geometry.solves,self.geometry.split_electric(weak))])

    def average(self, values):
        return np.sum(self.geometry.weight*values,axis=(2,3))/self.mesh.volumes

    def current_energy(self, current):
        j=self.geometry.sample_electric(current)
        return self.average(.5*self.inertia_coefficient*np.sum(j*j,axis=-1))

    def density_forcing(self, current, density_rate):
        """Weak -mu J n_t/(e² n²), using actual supplied cell-average n_t."""
        j=self.geometry.sample_electric(current)
        ratio=self.density_map.sample(density_rate)/self.density
        return weak_vector(self.geometry,-(self.inertia_coefficient*ratio)[...,None]*j)

    def density_energy_gradient(self, current):
        """Derivative of total K_J with respect to number-density cell means."""
        j=self.geometry.sample_electric(current)
        return self.density_map.adjoint(-.5*self.inertia_coefficient*
            np.sum(j*j,axis=-1)/self.density)

    def current_energy_rate(self, current, current_rate, density_rate):
        j=self.geometry.sample_electric(current);jt=self.geometry.sample_electric(current_rate)
        ratio=self.density_map.sample(density_rate)/self.density
        changing_current=self.average(self.inertia_coefficient*np.sum(j*jt,axis=-1))
        changing_density=self.average(-.5*self.inertia_coefficient*np.sum(j*j,axis=-1)*ratio)
        return dict(current=changing_current,density=changing_density,
                    total=changing_current+changing_density)

    def magnetic_energy_rate(self, field, magnetic_rate):
        """Full local B dot B_t, including harmonic bias cross work."""
        psi,g=self.geometry.unpack(magnetic_rate)
        if np.any(psi[[0,-1]]!=0):raise ValueError('Fixed axis/wall magnetic rate required')
        rate=self.geometry.magnetic_samples(psi,g)
        return self.average(np.sum(field.samples(self.geometry)*rate,axis=-1)/MU0)
