"""Ideal common-velocity MHD with separate advected electron and ion temperatures."""
import numpy as np
from .toroidal_mhd import ToroidalMHD
from .toroidal_state import ToroidalState
from .two_temperature import SpeciesThermodynamics
from .bulk_entropy import entropy_rate


class AdiabaticSpeciesMHD(ToroidalMHD):
    """Bounded ideal verification path, not yet the warm Hall/formation solver."""
    def __init__(self, mesh, initial, electron_pressure, ion_mass, parameters=None, *, density_scale=1e20):
        super().__init__(mesh, initial, parameters)
        if self.parameters.resistivity != 0:
            raise ValueError('Species heating partition is not implemented; ideal resistivity=0 required')
        self.thermal = SpeciesThermodynamics(mesh, ion_mass, self.parameters.gamma, density_scale)
        self.electron_entropy = self.thermal.from_pressure(self.state, electron_pressure)
        self.initial_species = self.thermal.diagnostics(self.state, self.electron_entropy)

    def species_rhs(self, state, entropy):
        self.thermal.validate(state, entropy)
        derivative, heat, statistics, mass_flux = self._rhs_with_fluxes(state)
        mapping = self.thermal.mapping
        flux = mapping.fluxes(*(v/self.thermal.ion_mass for v in mass_flux))
        du = entropy_rate(mapping, self.thermal.density(state), entropy, *flux, order=self.parameters.spatial_order)
        return derivative, du, heat, statistics

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f'Time step must be positive and <= fluid limit {limit:.6g} s')
        old, u = self.state, self.electron_entropy
        a, ua, qa, sa = self.species_rhs(old, u)
        predictor = ToroidalState(old.conserved+dt*a.conserved, old.psi+dt*a.psi, old.magnetic_quadrature)
        up = u+dt*ua
        b, ub, qb, sb = self.species_rhs(predictor, up)
        new = ToroidalState(old.conserved+.5*dt*(a.conserved+b.conserved), old.psi+.5*dt*(a.psi+b.psi), old.magnetic_quadrature)
        unew = u+.5*dt*(ua+ub)
        self.thermal.validate(new, unew)
        self.state, self.electron_entropy = new, unew
        self.time += dt
        self.steps += 1
        self.physical_ohmic_heat += .5*dt*(qa+qb)
        for key in self.limiter_totals:
            if key == 'minimum_pressure_blend_weight':
                self.limiter_totals[key] = min(self.limiter_totals[key], sa[key], sb[key])
            else:
                self.limiter_totals[key] += sa[key]+sb[key]
        return dt

    def diagnostics(self):
        result = super().diagnostics()
        species = self.thermal.diagnostics(self.state, self.electron_entropy)
        result.update(species)
        for name in ('electron_entropy_proxy_integral_j', 'electron_log_entropy_particles'):
            result[name+'_change'] = species[name]-self.initial_species[name]
        result['species_model'] = 'equal-gamma ideal common-velocity; no Hall, collisions or conduction'
        return result
