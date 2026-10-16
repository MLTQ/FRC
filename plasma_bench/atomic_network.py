"""Joint level/charge evolution with independent electron number and energy."""
import numpy as np
from scipy.constants import elementary_charge as EV, h, c
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from .adf07_ionization import coefficient


class AtomicNetwork:
    """Finite supplied atom; omitted channels and ground-collapse assumptions are explicit."""

    def __init__(self, ionization, recombination, levels=None, *, resolved_ionization=(),
                 ionization_enabled=True, recombination_enabled=True, radiation=True):
        self.ionization = ionization; self.recombination = recombination
        if ionization.nuclear_charge != recombination.nuclear_charge:
            raise ValueError('Matching complete charge ladders required')
        self.levels = dict(levels or {})
        self.ionization_enabled = ionization_enabled
        self.recombination_enabled = recombination_enabled
        self.radiation = radiation
        charges = []; energies = []; blocks = {}; ground = []; excited = []
        for charge in ionization.charge:
            level = self.levels.get(charge)
            if level is not None and (level.data['charge'] != charge or
                    level.data['atomic_number'] != ionization.nuclear_charge):
                raise ValueError('Resolved states must match their supplied element and charge')
            relative = np.zeros(1) if level is None else level.energy_ev
            start = len(charges); ground.append(start)
            blocks[charge] = np.arange(start, start+len(relative))
            charges.extend([charge]*len(relative))
            energies.extend(ionization.atomic_energy[charge]/EV+relative)
            excited.extend(range(start+1, start+len(relative)))
        if set(self.levels)-set(blocks):
            raise ValueError('Resolved charge outside complete ladder')
        self.charge = np.array(charges); self.energy_ev = np.array(energies)
        self.blocks = blocks; self.ground = np.array(ground); self.excited = np.array(excited, int)
        self.count = len(charges)
        self.electrons, self.thermal, self.photons, self.heating, self.constant = range(self.count, self.count+5)
        self.resolved_ionization = {}
        for charge, local_index, row in resolved_ionization:
            if (charge not in self.levels or not 0 <= local_index < len(blocks[charge]) or
                    row['charge'] != charge or row['final_charge'] != charge+1 or
                    row['final_metastable'] != 0):
                raise ValueError('Explicit matched initial state and unresolved adjacent destination required')
            origin = blocks[charge][local_index]
            if origin in self.resolved_ionization:
                raise ValueError('Multiple ionization tables for one origin would double count')
            self.resolved_ionization[origin] = row
        origins = sorted(set(self.ground[:-1]) | set(self.resolved_ionization))
        self.ion_origins = np.array(origins)
        self.ion_destinations = self.ground[self.charge[self.ion_origins]+1]
        self.ion_delta_ev = self.energy_ev[self.ion_destinations]-self.energy_ev[self.ion_origins]
        if np.any(self.ion_delta_ev <= 0):
            raise ValueError('Bound initial states and positive physical ionization energies required')
        self.missing_ionization = np.setdiff1d(self.excited, self.ion_origins)

    def initial(self, charge_fractions, temperature_ev):
        fractions = np.asarray(charge_fractions, float)
        if fractions.shape != self.ground.shape or np.any(fractions < 0) or not np.isfinite(fractions).all():
            raise ValueError('Nonnegative finite initial charge fractions required')
        state = np.zeros(self.count+4); state[self.ground] = fractions
        state[self.electrons] = self.charge@state[:self.count]
        state[self.thermal] = 1.5*state[self.electrons]*temperature_ev
        self.validate(state)
        return state

    def validate(self, state):
        if (state.shape != (self.count+4,) or not np.isfinite(state).all() or np.any(state < 0)):
            raise ValueError('Finite nonnegative populations and physical reservoirs required; no repair')
        electron = state[self.electrons]; energy = state[self.thermal]
        if (electron == 0 and energy != 0) or (electron > 0 and energy <= 0):
            raise ValueError('Electron energy must correspond to actual carriers, without a seed floor')

    def generator(self, state, density, heating_ev_per_nucleus_s=0.):
        self.validate(state)
        if not np.isfinite([density, heating_ev_per_nucleus_s]).all() or density <= 0 or heating_ev_per_nucleus_s < 0:
            raise ValueError('Positive nucleus density and nonnegative declared heat required')
        number = state[self.electrons]; matrix = np.zeros((self.count+5, self.count+5))
        # Spontaneous decay does not require a free-electron population.
        if self.radiation:
            for charge, level in self.levels.items():
                radiation = level.radiation; block = self.blocks[charge]
                lower, upper = block[radiation.lower], block[radiation.upper]
                np.add.at(matrix, (lower, upper), radiation.A)
                np.add.at(matrix, (upper, upper), -radiation.A)
                np.add.at(matrix[self.photons], upper, radiation.A*radiation.photon_energy/EV)
        if number == 0:
            if heating_ev_per_nucleus_s:
                raise ValueError('Thermal electron heating requires actual seed electrons')
            return matrix, dict(temperature_ev=None, electron_density_m3=0.,
                ionization_events_per_nucleus_s=np.zeros(len(self.ion_origins)),
                resolved_excited_ionization_events_per_nucleus_s=0.,
                recombination_events_per_nucleus_s=np.zeros(len(self.ground)-1),
                missing_excited_ionization_population=float(state[self.missing_ionization].sum()),
                missing_excited_parent_rr_population=float(state[self.excited[self.charge[self.excited]>0]].sum()),
                outside_ionization_fit_event_fraction=0., ionization_threshold_discrepancy_max_ev=0.,
                ionization_threshold_discrepancy_over_5pct_event_fraction=0.,
                unknown_rr_range_event_fraction=0., outside_known_rr_range_event_fraction=0.)
        temperature = state[self.thermal]/(1.5*number); ne = density*number
        ion_rates, inside = self.ionization.coefficients(temperature*EV)
        rr_rates, capture, rr_inside = self.recombination.coefficients(temperature*EV)
        for charge, level in self.levels.items():
            block = self.blocks[charge]; n = len(block)
            local = level.generator(temperature, ne, radiation=False)
            matrix[np.ix_(block, block)] += local[:n, :n]
            matrix[self.thermal, block] += local[n, :n]
            matrix[self.photons, block] += local[n+1, :n]
        rates = ne*ion_rates[self.charge[self.ion_origins]]
        valid = inside[self.charge[self.ion_origins]].copy()
        source_threshold = self.ionization.fit_energy[self.charge[self.ion_origins]]/EV
        for k, origin in enumerate(self.ion_origins):
            if origin in self.resolved_ionization:
                row = self.resolved_ionization[origin]
                rates[k] = ne*coefficient(row, temperature); valid[k] = True
                source_threshold[k] = row['effective_potential_cm1']*100*h*c/EV
        if not self.ionization_enabled:
            rates *= 0
        for origin, destination, rate, delta in zip(self.ion_origins, self.ion_destinations, rates, self.ion_delta_ev):
            matrix[destination, origin] += rate; matrix[origin, origin] -= rate
            matrix[self.electrons, origin] += rate
            matrix[self.thermal, origin] -= delta*rate
        rr_origins, rr_destinations = self.ground[1:], self.ground[:-1]
        rr_frequency = ne*rr_rates if self.recombination_enabled else np.zeros_like(rr_rates)
        for origin, destination, rate, loss in zip(rr_origins, rr_destinations, rr_frequency, capture/EV):
            matrix[destination, origin] += rate; matrix[origin, origin] -= rate
            matrix[self.electrons, origin] -= rate; matrix[self.thermal, origin] -= loss*rate
            matrix[self.photons, origin] += (self.energy_ev[origin]-self.energy_ev[destination]+loss)*rate
        matrix[self.thermal, self.constant] = heating_ev_per_nucleus_s
        matrix[self.heating, self.constant] = heating_ev_per_nucleus_s
        events = rates*state[self.ion_origins]; rr_events = rr_frequency*state[rr_origins]
        discrepancy = abs(source_threshold/self.ion_delta_ev-1)
        total = events.sum(); rr_total = rr_events.sum(); known = self.recombination.range_known
        diagnostics = dict(temperature_ev=float(temperature), electron_density_m3=float(ne),
            ionization_events_per_nucleus_s=events,
            resolved_excited_ionization_events_per_nucleus_s=float(events[self.ion_origins != self.ground[self.charge[self.ion_origins]]].sum()),
            recombination_events_per_nucleus_s=rr_events,
            missing_excited_ionization_population=float(state[self.missing_ionization].sum()),
            missing_excited_parent_rr_population=float(state[self.excited[self.charge[self.excited]>0]].sum()),
            outside_ionization_fit_event_fraction=float(events[~valid].sum()/total) if total else 0.,
            ionization_threshold_discrepancy_max_ev=float(np.max(abs(source_threshold-self.ion_delta_ev))),
            ionization_threshold_discrepancy_over_5pct_event_fraction=float(events[discrepancy>.05].sum()/total) if total else 0.,
            unknown_rr_range_event_fraction=float(rr_events[~known].sum()/rr_total) if rr_total else 0.,
            outside_known_rr_range_event_fraction=float(rr_events[known & ~rr_inside].sum()/rr_total) if rr_total else 0.)
        return matrix, diagnostics

    def step(self, state, dt, density, *, initial_heat=0., midpoint_heat=0.):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite timestep required')
        def advance(matrix, initial, duration):
            # The extra constant is an affine-source coordinate, not physical state.
            # Sparse action retains disconnected populations exactly; dense Padé
            # factorization introduced negative cross-block roundoff in pure decay.
            operator = csr_matrix(duration*matrix)
            return expm_multiply(operator, np.r_[initial, 1.],
                traceA=float(operator.diagonal().sum()))[:-1]
        first, _ = self.generator(state, density, initial_heat)
        midpoint = advance(first, state, .5*dt)
        second, _ = self.generator(midpoint, density, midpoint_heat)
        evolved = advance(second, state, dt)
        self.validate(evolved)
        # Validate all data intervals at the accepted endpoint as well as stages.
        self.generator(evolved, density, midpoint_heat)
        return evolved
