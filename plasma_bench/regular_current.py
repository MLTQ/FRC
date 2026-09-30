"""Physical current inertia and ideal/Hall work on the regular fluid state."""
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME
from scipy.sparse import block_diag
from scipy.sparse.linalg import LinearOperator
from .regular_field_map import RegularDensity
from .variable_compatible_current import VariableCompatibleCurrent,weak_vector


class RegularCurrent(VariableCompatibleCurrent):
    def __init__(self,mapping,number_density,ion_mass):
        if not np.isfinite(ion_mass) or ion_mass<=0:raise ValueError('Positive finite physical ion mass required')
        if not np.isclose(mapping.fluid.mass,ion_mass+ME,rtol=1e-14,atol=0.):raise ValueError('Same physical fluid and current masses required')
        self.mapping=mapping;self.geometry=mapping.geometry;self.mesh=self.geometry.mesh
        self.density_map=RegularDensity(mapping,number_density);self.density=self.density_map.samples
        self.ion_mass=ion_mass;self.electron_mass=ME;self.alpha=ME/(ion_mass+ME);self.beta=ion_mass/(ion_mass+ME)
        self.delta=self.beta-self.alpha;self.reduced_mass=ion_mass*self.alpha
        with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
            self.inertia_coefficient=self.reduced_mass/(self.density*E**2)
        if not np.isfinite(self.inertia_coefficient).all() or np.any(self.inertia_coefficient<=0):raise ValueError('Physical current inertia exceeds positive finite range')
        self.inertia_coefficient.setflags(write=False)
        self.mass=block_diag(self.geometry.mass,format='csc');self.hessian=self.geometry.hessian
        self.inertia=LinearOperator(self.mass.shape,matvec=self._inertia,rmatvec=self._inertia,dtype=float)

    def electromagnetic_pair(self,state,field,density_rate):
        """Instantaneous ideal/Hall pair with the explicit incomplete density exchange."""
        if not np.array_equal(state[...,0],self.density_map.point_values):raise ValueError('Current model must use this actual state density')
        geo=self.geometry;j=self.current(field.packed_increment(geo));physical_j=geo.sample_electric(j)
        v=self.mapping.velocity(state);b=field.samples(geo);force=np.cross(physical_j,b)
        ideal=weak_vector(geo,-np.cross(v,b));hall=weak_vector(geo,self.delta*force/(E*self.density[...,None]))
        changing=self.density_forcing(j,density_rate);rate=self.mapping.force_rate(state,force)
        exchange=float(-np.sum(self.density_energy_gradient(j)*density_rate))
        return dict(current=j,fluid_rate=rate,ideal_forcing=ideal,hall_forcing=hall,
                    density_forcing=changing,forcing=ideal+hall+changing,
                    unresolved_density_exchange=exchange)
