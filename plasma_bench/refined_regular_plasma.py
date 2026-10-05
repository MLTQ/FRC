"""Opt-in coupled assembly with full-residual-controlled induction refinement."""
from copy import copy
from .regular_plasma import RegularPlasma
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies
from .refined_kernel_induction import RefinedKernelInduction


class RefinedRegularPlasma(RegularPlasma):
    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model),
                         RefinedKernelInduction(model, **self.solver_options))
            if not self.cache_geometry:
                return operators
            self._operators = operators
        operators = tuple(copy(op) for op in self._operators)
        for op in operators:
            op.model = model
        return operators
