"""Original current-work expressions with expired quadrature temporaries released."""
import numpy as np
from scipy.constants import elementary_charge as E
from .regular_current import RegularCurrent
from .variable_compatible_current import weak_vector


class ReleasedRegularCurrent(RegularCurrent):
    def electromagnetic_pair(self,state,field,density_rate):
        """Instantaneous ideal/Hall pair with the explicit incomplete density exchange."""
        if not np.array_equal(state[...,0],self.density_map.point_values):raise ValueError('Current model must use this actual state density')
        geo=self.geometry;j=self.current(field.packed_increment(geo));physical_j=geo.sample_electric(j)
        v=self.mapping.velocity(state);b=field.samples(geo);force=np.cross(physical_j,b)
        ideal=weak_vector(geo,-np.cross(v,b));hall=weak_vector(geo,self.delta*force/(E*self.density[...,None]))
        del physical_j, v, b
        changing=self.density_forcing(j,density_rate);rate=self.mapping.force_rate(state,force)
        del force
        exchange=float(-np.sum(self.density_energy_gradient(j)*density_rate))
        return dict(current=j,fluid_rate=rate,ideal_forcing=ideal,hall_forcing=hall,
                    density_forcing=changing,forcing=ideal+hall+changing,
                    unresolved_density_exchange=exchange)
