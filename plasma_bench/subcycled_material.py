"""Explicit separation of source intervals and resolved material transport."""
import numpy as np
from .precise_reactions import PreciseReactiveMaterial


class SubcycledReactiveMaterial(PreciseReactiveMaterial):
    """Reaction half / many actual transport steps / reaction half."""
    def __init__(self, source, reference_density, *, transport_fraction=.5,
                 maximum_transport_steps=100000, **controls):
        super().__init__(source, reference_density, **controls)
        if not np.isfinite(transport_fraction) or not 0 < transport_fraction <= 1:
            raise ValueError('Transport fraction must be positive and no larger than one')
        if not isinstance(maximum_transport_steps, int) or maximum_transport_steps < 1:
            raise ValueError('Positive explicit transport-step budget required')
        self.transport_fraction = float(transport_fraction)
        self.maximum_transport_steps = maximum_transport_steps
        self.last_transport = []

    def step(self, model, state, reservoirs, dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite source interval required')
        self.last_transport = []
        model, state, reservoirs = self.source_step(model, state, reservoirs, dt/2)
        elapsed = 0.
        while elapsed < dt:
            if len(self.last_transport) >= self.maximum_transport_steps:
                raise ValueError('Declared transport-step budget; source interval not complete')
            bound = model.timestep(state)
            h = min(self.transport_fraction*bound, dt-elapsed)
            if not np.isfinite(h) or h <= 0 or elapsed+h == elapsed:
                raise ValueError('Transport step cannot advance represented time')
            state = model.step(state, h)
            elapsed = min(elapsed+h, dt)
            budget = self.totals(model, state, reservoirs)
            self.last_transport.append(dict(time_within_transport_s=elapsed, dt_s=h,
                admissible_dt_s=bound, energy_j_m2=budget['full_energy_j_m2'],
                electric_energy_j_m2=budget['electric_energy_j_m2'],
                minimum_thermal_j_m2=min(float(s.thermal.min()*m.number_scale)
                                        for m, s in zip(model.models, state.blocks))))
        return self.source_step(model, state, reservoirs, dt/2)
