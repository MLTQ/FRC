"""Prescribed physical resistivity and paired heating in the moving wall plasma."""
from copy import copy
from dataclasses import replace
import numpy as np
from .current_transport import CurrentTransport
from .regular_current import RegularCurrent
from .regular_species import RegularSpecies
from .split_inverse_wall_induction import SplitInverseWallInduction
from .split_inverse_wall_plasma import SplitInverseWallPlasma
from .variable_compatible_current import weak_vector


class ResistiveTransport(CurrentTransport):
    """Keep all current transport, adding eta*J and its local thermal work."""
    def __init__(self, model, resistivity, electron_heat_fraction):
        super().__init__(model)
        self.resistivity = resistivity
        self.fractions = np.array([electron_heat_fraction, 1-electron_heat_fraction])

    def pair(self, state, current, density_rate):
        result = super().pair(state, current, density_rate)
        model = self.model
        j = model.geometry.sample_electric(current)
        heat = self.resistivity*np.sum(j*j, axis=-1)
        # V*h*Cdot = projection of Q: the same physical Joule work as eta*J.
        nodal_power = model.mapping.adjoint(heat)/model.mapping.fluid.volume
        h = model.mapping.fluid.fields(state)['specific_internal']
        result['fluid_rate'][..., 4:] += nodal_power[..., None]*self.fractions/h
        result['forcing'] += weak_vector(model.geometry, self.resistivity*j)
        if not all(np.isfinite(result[k]).all() for k in ('fluid_rate', 'forcing')):
            raise ValueError('Nonfinite resistive force or thermal source')
        return result


class ResistiveWallPlasma(SplitInverseWallPlasma):
    """Full moving equations with explicitly prescribed scalar eta and heat split."""
    def __init__(self, mapping, ion_mass, *, resistivity, electron_heat_fraction, **options):
        if not np.isfinite(resistivity) or resistivity < 0:
            raise ValueError('Finite nonnegative resistivity in ohm meter required')
        if not np.isfinite(electron_heat_fraction) or not 0 <= electron_heat_fraction <= 1:
            raise ValueError('Explicit electron heating fraction in [0,1] required')
        super().__init__(mapping, ion_mass, **options)
        self.resistivity = float(resistivity)
        self.electron_heat_fraction = float(electron_heat_fraction)

    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (ResistiveTransport(model, self.resistivity, self.electron_heat_fraction),
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
        j = self.geometry.sample_electric(model.current(field.packed_increment(self.geometry)))
        q = self.resistivity*np.sum(j*j, axis=-1)
        diagnostics = dict(result.diagnostics,
            ohmic_power=float(np.sum(self.geometry.weight*q)),
            resistivity_ohm_m=self.resistivity,
            electron_heat_fraction=self.electron_heat_fraction)
        return replace(result, diagnostics=diagnostics)
