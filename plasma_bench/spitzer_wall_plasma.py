"""Weak-field Z=1 electrical resistance, Joule heat and temperature exchange."""
from copy import copy
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0
from .coulomb_exchange_wall_plasma import CoulombExchangeWallPlasma, equilibration_rate
from .current_transport import CurrentTransport
from .regular_current import RegularCurrent
from .regular_species import RegularSpecies
from .split_inverse_wall_induction import SplitInverseWallInduction
from .variable_compatible_current import weak_vector


def spitzer_samples(mapping, state, ion_mass):
    """Actual physical temperature and conditional low-field coefficients."""
    fluid = mapping.fluid
    fluid.fields(state)
    n = mapping.scalar(state[..., 0])
    if np.any(n <= 0):
        raise ValueError('Positive reconstructed plasma density required')
    entropy = np.stack([mapping.scalar(state[..., k])/n for k in (4, 5)], axis=-1)
    temperature = np.exp(entropy+(fluid.thermal.gammas-1)*np.log(n[..., None]/fluid.thermal.n0))
    te, ti = np.moveaxis(temperature, -1, 0)
    _, log = equilibration_rate(n, te, ti, ion_mass)
    coupling = E**4/(4*np.pi*EPS0)**2
    frequency = 4*np.sqrt(2*np.pi)*n*coupling*log/(3*np.sqrt(ME)*(E*te)**1.5)
    eta = .51*ME*frequency/(n*E**2)
    if not np.isfinite(eta).all() or np.any(eta <= 0):
        raise ValueError('Nonpositive or nonfinite Spitzer resistivity')
    return n, temperature, eta, frequency, log


class SpitzerTransport(CurrentTransport):
    def pair(self, state, current, density_rate):
        result = super().pair(state, current, density_rate)
        model = self.model
        _, _, eta, _, _ = spitzer_samples(model.mapping, state, model.ion_mass)
        j = model.geometry.sample_electric(current)
        heat = eta*np.sum(j*j, axis=-1)
        nodal_power = model.mapping.adjoint(heat)/model.mapping.fluid.volume
        h = model.mapping.fluid.fields(state)['specific_internal']
        result['fluid_rate'][..., 4] += nodal_power/h[..., 0]
        result['forcing'] += weak_vector(model.geometry, eta[..., None]*j)
        if not all(np.isfinite(result[k]).all() for k in ('fluid_rate', 'forcing')):
            raise ValueError('Nonfinite Spitzer force or electron heat source')
        return result


class SpitzerWallPlasma(CoulombExchangeWallPlasma):
    """Conditional weak-field closure; no complete Braginskii transport claim."""
    def __init__(self, mapping, ion_mass, **options):
        # Parent supplies thermal exchange only; variable eta replaces its operator.
        super().__init__(mapping, ion_mass, resistivity=0., electron_heat_fraction=1., **options)

    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (SpitzerTransport(model), RegularSpecies(model),
                SplitInverseWallInduction(model, **self.solver_options))
            if not self.cache_geometry:
                return operators
            self._operators = operators
        operators = tuple(copy(operator) for operator in self._operators)
        for operator in operators:
            operator.model = model
        return operators

    def derivative(self, state, field):
        result = super().derivative(state, field)
        model = RegularCurrent(self.mapping, state[..., 0], self.ion_mass)
        packed = field.packed_increment(self.geometry)
        j = self.geometry.sample_electric(model.current(packed))
        b = field.samples(self.geometry)
        n, temperature, eta, frequency, log = spitzer_samples(self.mapping, state, self.ion_mass)
        hall = E*np.linalg.norm(b, axis=-1)/(ME*frequency)
        drift = np.linalg.norm(j, axis=-1)/(n*E)
        drift_ratio = drift[..., None]/np.sqrt(E*temperature/np.array([ME, self.ion_mass]))
        diagnostics = dict(result.diagnostics)
        diagnostics.pop('resistivity_ohm_m')
        diagnostics.update(ohmic_power=float(np.sum(self.geometry.weight*eta*np.sum(j*j, axis=-1))),
            resistivity_range_ohm_m=[float(eta.min()), float(eta.max())],
            electron_collision_rate_range_s=[float(frequency.min()), float(frequency.max())],
            spitzer_coulomb_log_range=[float(log.min()), float(log.max())],
            maximum_electron_hall_parameter=float(hall.max()),
            maximum_drift_over_species_thermal_speed=np.max(drift_ratio.reshape(-1, 2), axis=0).tolist(),
            weak_field_review=bool(np.any(hall >= .1)),
            relative_drift_review=bool(np.any(drift_ratio >= .1)),
            spitzer_scope='Z1 weak-field local Maxwellian friction and electron Joule heat; thermal exchange included, thermal forces/conduction/viscosity/neutrals absent. Review thresholds are not validity proofs.')
        return replace(result, diagnostics=diagnostics)
