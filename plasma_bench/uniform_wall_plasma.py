"""Unchanged wall equations with physical-current-inertia preconditioning."""
from copy import copy
from .strong_wall_induction import StrongWallInduction
from .strong_wall_plasma import StrongWallPlasma
from .uniform_wall_metric import UniformWallMetric
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies


class UniformWallInduction(StrongWallInduction):
    def __init__(self, model, **options):
        super().__init__(model, **options)
        self.factor = UniformWallMetric(model, self.lift)


class UniformWallPlasma(StrongWallPlasma):
    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model), UniformWallInduction(model, **self.solver_options))
            if not self.cache_geometry: return operators
            self._operators = operators
        operators = tuple(copy(operator) for operator in self._operators)
        for operator in operators: operator.model = model
        return operators


def preconditioned(system):
    result = UniformWallPlasma(system.mapping, system.ion_mass,
        cache_geometry=system.cache_geometry, **system.solver_options)
    result.original_system = system.original_system
    return result
