"""Opt-in sixth-order numerical dissipation with explicitly paired heat."""
from copy import copy
from dataclasses import replace
import numpy as np
from scipy.constants import mu_0 as MU0
from .resistive_wall_plasma import ResistiveTransport
from .regular_current import RegularCurrent
from .regular_species import RegularSpecies
from .split_inverse_wall_induction import SplitInverseWallInduction
from .split_inverse_wall_plasma import SplitInverseWallPlasma


def current_laplacian(model, current):
    """Positive, electric-mass self-adjoint K = -mu0 Ampere Faraday."""
    return -MU0*model.current(model.geometry.faraday@current)


class HyperdiffusiveTransport(ResistiveTransport):
    def __init__(self, model, eta6, electron_heat_fraction):
        super().__init__(model, 0., electron_heat_fraction)
        self.eta6 = eta6

    def pair(self, state, current, density_rate):
        result = super().pair(state, current, density_rate)
        model = self.model
        kj = current_laplacian(model, current)
        result['forcing'] += self.eta6*(model.mass@current_laplacian(model, kj))
        q = self.eta6*np.sum(model.geometry.sample_electric(kj)**2, axis=-1)
        nodal_power = model.mapping.adjoint(q)/model.mapping.fluid.volume
        h = model.mapping.fluid.fields(state)['specific_internal']
        result['fluid_rate'][..., 4:] += nodal_power[..., None]*self.fractions/h
        return result


class HyperdiffusiveWallPlasma(SplitInverseWallPlasma):
    def __init__(self, mapping, ion_mass, *, eta6, electron_heat_fraction, **options):
        if not np.isfinite(eta6) or eta6 < 0:
            raise ValueError('Finite nonnegative numerical eta6 in ohm meter^5 required')
        if not np.isfinite(electron_heat_fraction) or not 0 <= electron_heat_fraction <= 1:
            raise ValueError('Explicit numerical electron heating fraction required')
        super().__init__(mapping, ion_mass, **options)
        self.eta6 = float(eta6)
        self.electron_heat_fraction = float(electron_heat_fraction)

    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (HyperdiffusiveTransport(model, self.eta6, self.electron_heat_fraction),
                RegularSpecies(model), SplitInverseWallInduction(model, **self.solver_options))
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
        current = model.current(field.packed_increment(self.geometry))
        kj = self.geometry.sample_electric(current_laplacian(model, current))
        power = float(self.eta6*np.sum(self.geometry.weight*np.sum(kj*kj, axis=-1)))
        return replace(result, diagnostics=dict(result.diagnostics,
            numerical_hyper_power=power, numerical_eta6_ohm_m5=self.eta6,
            numerical_electron_heat_fraction=self.electron_heat_fraction))
