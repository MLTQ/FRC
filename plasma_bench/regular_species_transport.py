"""Atomic fixed-current RK4 verification of regular species entropy transport."""
import numpy as np
from scipy.constants import elementary_charge as E
from .regular_species import RegularSpecies


class FrozenSpeciesTransport:
    """Density and current are prescribed, not advanced or called equilibrium."""
    def __init__(self, model, current, species=0):
        if species not in (0, 1):
            raise ValueError('Electron (0) or ion (1) species required')
        self.pair = RegularSpecies(model)
        self.model = model
        self.fluid = model.mapping.fluid
        self.n = model.density_map.point_values
        self.current = np.asarray(current, float).copy()
        if self.current.shape != (model.mass.shape[0],) or not np.isfinite(self.current).all():
            raise ValueError('Finite complete prescribed current required')
        self.current.setflags(write=False)
        self.fraction = (-model.beta, model.alpha)[species]
        self.gamma = self.fluid.thermal.gammas[species]

    def validate(self, entropy):
        entropy = self.model.mapping.validate(entropy)
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            h = E*np.exp(entropy+(self.gamma-1)*np.log(self.n/self.fluid.thermal.n0))/(self.gamma-1)
        if not np.isfinite(h).all() or np.any(h <= 0):
            raise ValueError('Positive finite species internal energy required; no repair')
        return entropy

    def rate(self, entropy):
        entropy = self.validate(entropy)
        return -self.fraction/E*self.pair.transport(self.current, entropy)/(self.fluid.volume*self.n)

    def step(self, entropy, dt):
        s = self.validate(entropy)
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite transport step required')
        k1 = self.rate(s)
        k2 = self.rate(s+.5*dt*k1)
        k3 = self.rate(s+.5*dt*k2)
        k4 = self.rate(s+dt*k3)
        result = s+dt*(k1+2*k2+2*k3+k4)/6
        self.validate(result)
        return result
