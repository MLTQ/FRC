"""Reactive atomic populations with separately evolved species momentum and energy."""
import numpy as np
from scipy.constants import elementary_charge as EV, electron_mass as ME
from scipy.sparse.linalg import expm_multiply


class AtomicMoments:
    """Homogeneous small-drift reaction moments; no spatial/field coupling yet."""
    def __init__(self, network, neutral_mass):
        self.network = network; self.neutral_mass = float(neutral_mass)
        self.mass = np.r_[self.neutral_mass-network.charge*ME, ME]
        if not np.isfinite(self.mass).all() or np.any(self.mass <= 0):
            raise ValueError('Positive physical neutral and charge-state masses required')
        self.speed_scale = np.sqrt(EV/self.neutral_mass)
        self.momentum_scale = self.neutral_mass*self.speed_scale
        self.base_size = network.count+4
        self.momentum = np.arange(self.base_size, self.base_size+3*(network.count+1)).reshape(-1, 3)
        self.heavy_energy = np.arange(self.momentum[-1, -1]+1, self.momentum[-1, -1]+1+network.count)
        self.energy = np.r_[self.heavy_energy, network.thermal]
        self.size = int(self.heavy_energy[-1]+1)

    def initial(self, atomic_state, heavy_temperature_ev, velocity):
        n = self.network; n.validate(atomic_state)
        velocity = np.asarray(velocity, float)
        temperature = np.broadcast_to(heavy_temperature_ev, (n.count,))
        if velocity.shape != (n.count+1, 3) or not np.isfinite(velocity).all() or not np.isfinite(temperature).all() or np.any(temperature < 0):
            raise ValueError('Finite species velocities and nonnegative heavy temperatures required')
        state = np.zeros(self.size); state[:self.base_size] = atomic_state
        numbers = np.r_[atomic_state[:n.count], atomic_state[n.electrons]]
        state[self.momentum] = (self.mass[:, None]/self.momentum_scale*velocity)*numbers[:, None]
        state[self.heavy_energy] = 1.5*numbers[:-1]*temperature
        state[self.energy] += (.5*(self.mass/EV)*np.sum(velocity**2, axis=1))*numbers
        self.fields(state)
        return state

    def fields(self, state):
        n = self.network
        if state.shape != (self.size,) or not np.isfinite(state).all():
            raise ValueError('Finite complete moment state required')
        numbers = np.r_[state[:n.count], state[n.electrons]]
        if np.any(numbers < 0):raise ValueError('Negative population; no clipping')
        empty = numbers == 0; momentum = state[self.momentum]; energy = state[self.energy]
        if np.any(momentum[empty] != 0) or np.any(energy[empty] != 0):
            raise ValueError('Empty species must have exactly zero momentum and energy')
        velocity = np.zeros_like(momentum)
        # Divide the dimensionless moments first. Multiplying a tiny but
        # representable population by a kilogram mass underflowed to zero.
        velocity[~empty] = (momentum[~empty]/numbers[~empty, None])*(self.momentum_scale/self.mass[~empty, None])
        kinetic = (.5*(self.mass/EV)*np.sum(velocity**2, axis=1))*numbers
        thermal = energy-kinetic
        if not np.isfinite(thermal).all() or np.any(thermal < 0):
            raise ValueError('Nonfinite or negative species thermal energy; no energy repair')
        atomic = state[:self.base_size].copy(); atomic[n.thermal] = thermal[-1]; n.validate(atomic)
        temperature = np.divide(thermal, 1.5*numbers, out=np.zeros_like(thermal), where=~empty)
        return dict(numbers=numbers, velocity=velocity, kinetic_ev=kinetic,
                    thermal_ev=thermal, temperature_ev=temperature, atomic_state=atomic)

    def generator(self, state, density, heat=0.):
        n = self.network; f = self.fields(state)
        base, diagnostic = n.generator(f['atomic_state'], density, heat)
        matrix = np.zeros((self.size+1, self.size+1))
        matrix[:self.base_size, :self.base_size] = base[:-1, :-1]
        matrix[:self.base_size, -1] = base[:-1, -1]
        ve = f['velocity'][-1]; ve2 = float(ve@ve)
        captures = (n.recombination.coefficients(f['temperature_ev'][-1]*EV)[1]/EV
                    if f['numbers'][-1] else np.zeros(n.ionization.nuclear_charge))
        transitions = base[:n.count, :n.count].copy(); np.fill_diagonal(transitions, 0.)
        if np.any(transitions < 0):raise ValueError('Atomic transition rates must be nonnegative')
        for destination, origin in np.argwhere(transitions > 0):
            rate = transitions[destination, origin]
            before, after = self.mass[[origin, destination]]
            change = n.charge[destination]-n.charge[origin]
            po, pd = self.momentum[[origin, destination]]
            eo, ed = self.heavy_energy[[origin, destination]]
            matrix[np.ix_(po, po)] -= rate*np.eye(3)
            matrix[eo, eo] -= rate
            if change in (0, 1):
                ratio = after/before
                matrix[np.ix_(pd, po)] += rate*ratio*np.eye(3)
                matrix[ed, eo] += rate*ratio
                if change == 1:
                    matrix[np.ix_(self.momentum[-1], po)] += rate*ME/before*np.eye(3)
                    matrix[n.thermal, eo] += rate*ME/before
            elif change == -1:
                capture = captures[n.charge[origin]-1]
                electron_p = rate*ME*ve/self.momentum_scale
                matrix[np.ix_(pd, po)] += rate*np.eye(3)
                matrix[pd, origin] += electron_p
                matrix[self.momentum[-1], origin] -= electron_p
                # Captured pair forms a heavy particle at its centre-of-mass
                # velocity. Relative kinetic energy joins the emitted binding energy.
                matrix[ed, eo] += rate*before/after
                cross = rate*ME/after*self.momentum_scale/EV*ve
                matrix[ed, po] += cross
                matrix[ed, origin] += rate*(.5*ME**2/after*ve2/EV+ME/after*capture)
                matrix[n.thermal, origin] -= rate*.5*ME*ve2/EV
                matrix[n.photons, eo] += rate*ME/after
                matrix[n.photons, po] -= cross
                matrix[n.photons, origin] += rate*(.5*before*ME/after*ve2/EV-ME/after*capture)
            else:
                raise ValueError('Only adjacent charge changes have specified moment closure')
        occupied = f['numbers'][:-1] > 0
        if f['numbers'][-1]:
            electron_speed = np.sqrt(2*EV*f['temperature_ev'][-1]/ME)
            relative = np.linalg.norm(f['velocity'][:-1][occupied]-ve, axis=1)/electron_speed
            thermal_ratio = np.sqrt(ME*f['temperature_ev'][:-1][occupied]/
                                    (self.mass[:-1][occupied]*f['temperature_ev'][-1]))
            diagnostic.update(maximum_drift_over_electron_thermal=float(relative.max(initial=0.)),
                              maximum_heavy_thermal_over_electron_thermal=float(thermal_ratio.max(initial=0.)))
        else:
            diagnostic.update(maximum_drift_over_electron_thermal=None,
                              maximum_heavy_thermal_over_electron_thermal=None)
        return matrix, diagnostic

    def step(self, state, dt, density, *, initial_heat=0., midpoint_heat=0.):
        if not np.isfinite(dt) or dt <= 0:raise ValueError('Positive finite step required')
        def advance(matrix, initial, duration):
            operator = duration*matrix
            return expm_multiply(operator, np.r_[initial, 1.], traceA=float(operator.diagonal().sum()))[:-1]
        first, _ = self.generator(state, density, initial_heat)
        midpoint = advance(first, state, .5*dt)
        second, _ = self.generator(midpoint, density, midpoint_heat)
        result = advance(second, state, dt)
        self.generator(result, density, midpoint_heat)
        return result

    def observe(self, state, density):
        n = self.network; f = self.fields(state); numbers = f['numbers']
        velocity = f['velocity']; charge_velocity = np.sum(n.charge[:, None]*numbers[:-1, None]*velocity[:-1], axis=0)
        return dict(number=float(numbers[:-1].sum()), electrons=float(numbers[-1]),
            charge_balance=float(numbers[-1]-n.charge@numbers[:-1]),
            mass_over_neutral=float(self.mass@numbers/self.neutral_mass),
            momentum_scaled=state[self.momentum].sum(axis=0),
            energy_budget_ev=float(n.energy_ev@numbers[:-1]+state[self.energy].sum()+state[n.photons]-state[n.heating]),
            current_a_m2=EV*density*(charge_velocity-numbers[-1]*velocity[-1]),
            electron_temperature_ev=float(f['temperature_ev'][-1]),
            heavy_thermal_ev=float(f['thermal_ev'][:-1].sum()),
            electron_thermal_ev=float(f['thermal_ev'][-1]),
            photon_energy_ev=float(state[n.photons]))
