"""Opt-in complete wall plasma with unchanged strict dual-residual bounds."""
from copy import copy
from .strong_wall_plasma import StrongWallPlasma
from .balanced_wall_induction import BalancedWallInduction
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies


class BalancedWallPlasma(StrongWallPlasma):
    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model), BalancedWallInduction(model, **self.solver_options))
            if not self.cache_geometry: return operators
            self._operators = operators
        operators = tuple(copy(operator) for operator in self._operators)
        for operator in operators: operator.model = model
        return operators


def balanced(system):
    result = BalancedWallPlasma(system.mapping, system.ion_mass,
        cache_geometry=system.cache_geometry, **system.solver_options)
    result.original_system = system.original_system
    return result
