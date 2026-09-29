"""Thermodynamic state for equal-gamma adiabatic ions and massless electrons."""
import numpy as np
from scipy.constants import elementary_charge as E
from .dual_mass import DualMass
from .toroidal_state import primitive
from .mhd_state import NonphysicalState


class SpeciesThermodynamics:
    def __init__(self, mesh, ion_mass, gamma=5/3, density_scale=1e20):
        if not all(np.isfinite(v) and v > 0 for v in (ion_mass, density_scale)) or not np.isfinite(gamma) or gamma <= 1:
            raise ValueError('Positive finite mass/density scale and gamma > 1 required')
        self.mapping, self.ion_mass = DualMass(mesh), ion_mass
        self.gamma, self.density_scale = gamma, density_scale

    def density(self, state):
        return self.mapping.density(state.conserved[..., 0])/self.ion_mass

    def pressure(self, state, entropy):
        n = self.density(state)
        u = self.mapping.geometry.scalar(entropy, 'electron entropy proxy')
        if np.any(n <= 0) or np.any(u <= 0):
            raise NonphysicalState('Nonpositive electron density/entropy; no floor applied')
        return (n/self.density_scale)**(self.gamma-1)*u

    def from_pressure(self, state, pressure):
        p = self.mapping.geometry.scalar(pressure, 'electron pressure')
        u = p/(self.density(state)/self.density_scale)**(self.gamma-1)
        self.validate(state, u)
        return u.copy()

    def validate(self, state, entropy):
        q = primitive(self.mapping.mesh, state, self.gamma)
        pe = self.pressure(state, entropy)
        pi = q[..., 3]-self.mapping.primal_average(pe)
        if not np.isfinite(pe).all() or not np.isfinite(pi).all() or np.any(pi <= 0):
            raise NonphysicalState('Nonpositive ion/electron thermal partition; no total-energy repair applied')
        return pe, pi

    def pressure_rate(self, state, entropy, mass_rate, entropy_rate):
        n = self.density(state)
        dn = self.mapping.density(mass_rate)/self.ion_mass
        pe = self.pressure(state, entropy)
        return (n/self.density_scale)**(self.gamma-1)*entropy_rate+(self.gamma-1)*pe*dn/n

    def diagnostics(self, state, entropy):
        pe, pi = self.validate(state, entropy)
        n, volume = self.density(state), self.mapping.volumes
        return dict(electron_internal_energy_j=float(np.sum(volume*pe)/(self.gamma-1)),
                    ion_internal_energy_j=float(np.sum(self.mapping.mesh.volumes*pi)/(self.gamma-1)),
                    electron_entropy_proxy_integral_j=float(np.sum(volume*entropy)),
                    electron_log_entropy_particles=float(np.sum(volume*n*np.log(entropy/(n*E)))),
                    min_electron_pressure_pa=float(pe.min()), min_ion_pressure_pa=float(pi.min()),
                    min_electron_temperature_ev=float(np.min(pe/(n*E))),
                    min_ion_temperature_ev=float(np.min(pi/(state.conserved[..., 0]/self.ion_mass*E))))
