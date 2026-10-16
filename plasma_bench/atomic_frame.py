"""Fixed Galilean coordinates for the homogeneous reactive moment source."""
import numpy as np
from scipy.constants import elementary_charge as EV
from .atomic_moments import AtomicMoments


class FixedFrameAtomicMoments(AtomicMoments):
    """Evolve relative moments while accepting and reporting laboratory motion.

    ``fields`` retains the inherited integration-frame contract. Use
    ``laboratory_state`` for extensive lab moments or ``observe`` for lab budgets.
    The reference velocity stays fixed throughout a trajectory.
    """
    def __init__(self, network, neutral_mass, reference_velocity):
        super().__init__(network, neutral_mass)
        frame = np.asarray(reference_velocity, dtype=float).copy()
        if frame.shape != (3,) or not np.isfinite(frame).all():
            raise ValueError('A finite three-component fixed reference velocity is required')
        frame.setflags(write=False)
        self.reference_velocity = frame

    def initial(self, atomic_state, heavy_temperature_ev, velocity):
        return super().initial(atomic_state, heavy_temperature_ev,
                               np.asarray(velocity, dtype=float)-self.reference_velocity)

    def laboratory_state(self, state):
        """Return physical lab extensive moments; do not use these for stepping."""
        f = self.fields(state); u = self.reference_velocity
        result = state.copy()
        result[self.momentum] += (self.mass[:, None]/self.momentum_scale*u)*f['numbers'][:, None]
        result[self.energy] += (self.momentum_scale/EV)*(state[self.momentum]@u)
        result[self.energy] += (.5*(self.mass/EV)*(u@u))*f['numbers']
        return result

    def observe(self, state, density):
        row = super().observe(state, density)
        n = self.network; u = self.reference_velocity
        numbers = np.r_[state[:n.count], state[n.electrons]]
        total_mass = self.mass@numbers
        relative_momentum = state[self.momentum].sum(axis=0)
        row['momentum_scaled'] = relative_momentum+total_mass/self.momentum_scale*u
        row['energy_budget_ev'] += (self.momentum_scale/EV)*(relative_momentum@u)+.5*total_mass/EV*(u@u)
        # Keep any actual charge residual; do not project onto quasineutrality.
        row['current_a_m2'] -= EV*density*row['charge_balance']*u
        return row
