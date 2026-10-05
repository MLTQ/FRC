"""Bounded plasma candidate with admissible normal velocity and toroidal B."""
from copy import copy
from dataclasses import replace
import numpy as np
from .regular_coordinate_fluid import RegularCoordinateFluid
from .regular_plasma import RegularPlasma
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies
from .strong_wall_induction import StrongWallInduction


class StrongWallFluid(RegularCoordinateFluid):
    def fields(self, state):
        result = super().fields(state)
        if np.any(np.asarray(state)[-1, -1, :, 1] != 0):
            raise ValueError('Inactive zero wall-normal velocity slots required')
        return result

    def derivative(self, state):
        rate, _ = super().derivative(state); rate[-1, -1, :, 1] = 0.
        return rate, dict(physical_energy_rate=np.sum(self.fields(state)['energy_gradient']*rate, axis=-1))


class StrongWallPlasma(RegularPlasma):
    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model), StrongWallInduction(model, **self.solver_options))
            if not self.cache_geometry: return operators
            self._operators = operators
        operators = tuple(copy(op) for op in self._operators)
        for op in operators: op.model = model
        return operators

    def derivative(self, state, field):
        result = super().derivative(state, field)
        fluid = result.fluid.copy(); reaction = -fluid[-1, -1, :, 1].copy(); fluid[-1, -1, :, 1] = 0.
        diagnostics = dict(result.diagnostics, wall_normal_reaction_max=float(np.max(abs(reaction))))
        return replace(result, fluid=fluid, diagnostics=diagnostics)
