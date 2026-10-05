"""Unforced bounded moving plasma assembled on one actual regular fluid state."""
from dataclasses import dataclass
from copy import copy
import numpy as np
from .regular_current import RegularCurrent
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies
from .kernel_induction import KernelInduction


@dataclass(frozen=True)
class RegularPlasmaRate:
    fluid: np.ndarray
    magnetic: np.ndarray
    current: np.ndarray
    electric: np.ndarray
    residual: float
    diagnostics: dict


class RegularPlasma:
    """Filled periodic cylinder; physical masses, no sources or exterior model."""
    def __init__(self, mapping, ion_mass, *, tolerance=1e-10, max_iterations=1000, cache_geometry=True):
        self.mapping = mapping
        self.fluid = mapping.fluid
        self.geometry = mapping.geometry
        self.ion_mass = ion_mass
        self.solver_options = dict(tolerance=tolerance, max_iterations=max_iterations)
        self.cache_geometry = cache_geometry
        self._operators = None

    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model),
                         KernelInduction(model, **self.solver_options))
            if not self.cache_geometry:
                return operators
            self._operators = operators
        # Only geometry matrices/factors are shared. No density-dependent
        # coefficient, rate, inertia action or previous solution is reused.
        operators = tuple(copy(op) for op in self._operators)
        for op in operators:
            op.model = model
        return operators

    def derivative(self, state, field):
        state = np.asarray(state, float)
        fields = self.fluid.fields(state)
        field.packed_increment(self.geometry)
        model = RegularCurrent(self.mapping, state[..., 0], self.ion_mass)
        bulk, _ = self.fluid.derivative(state)
        nt = bulk[..., 0]
        electromagnetic = model.electromagnetic_pair(state, field, nt)
        current = electromagnetic['current']
        transport_operator, species_operator, induction = self.operators(model)
        transport = transport_operator.pair(state, current, nt)
        species = species_operator.pair(state, current)
        fluid = bulk + electromagnetic['fluid_rate'] + transport['fluid_rate'] + species['fluid_rate']
        forcing = electromagnetic['forcing'] + transport['forcing'] + species['reduced_forcing']
        solved = induction.solve(forcing)
        # Compatible species gradients are exact curl-kernel fields. Keep their
        # analytic split in evolution; restore the full E for original Ohm work.
        electric = solved.electric + species['potential_electric']
        full_forcing = electromagnetic['forcing'] + transport['forcing'] + species['forcing']
        inertia = model.inertia @ solved.current
        terms = abs(model.mass @ electric) + abs(inertia) + abs(full_forcing)
        backward = float(np.linalg.norm(model.mass @ electric - inertia - full_forcing)
                         / max(np.linalg.norm(terms), np.finfo(float).tiny))
        residual = max(solved.residual, backward)
        if not np.isfinite(residual) or residual > self.solver_options['tolerance']:
            raise RuntimeError('Full species-restored Ohm equation failed; no input modified')
        if not all(np.isfinite(a).all() for a in (fluid, electric, solved.magnetic, solved.current)):
            raise ValueError('Nonfinite coupled rate; no input modified')
        grid_fluid = float(np.sum(self.fluid.volume[..., None]*fields['energy_gradient']*fluid))
        current_power = model.current_energy_rate(current, solved.current, nt)['total']
        magnetic_power = model.magnetic_energy_rate(field, solved.magnetic)
        powers = dict(fluid=grid_fluid, current=float(np.sum(model.mesh.volumes*current_power)),
                      magnetic=float(np.sum(model.mesh.volumes*magnetic_power)))
        diagnostics = dict(powers=powers, total_power=sum(powers.values()),
                           full_species_ohm_backward=backward, iterations=solved.iterations,
                           induction_residuals=solved.residual_components)
        return RegularPlasmaRate(fluid, solved.magnetic, solved.current, electric, residual, diagnostics)
