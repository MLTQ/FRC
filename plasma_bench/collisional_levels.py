"""Resolved electron-impact excitation with electron and photon energy exchange."""
import numpy as np
from scipy.constants import elementary_charge as E, Boltzmann as KB, electron_mass as ME, physical_constants
from scipy.linalg import expm
from .adf04_levels import read_levels
from .line_radiation import RadiativeLevels


class CollisionalLevels:
    """Finite specified level set and Maxwellian bath; not a complete charge model."""
    def __init__(self, path, *, radiative_path=None, exclude_processes=()):
        self.data = read_levels(path); rows = self.data['transitions']
        supplied = {r['code'] for r in self.data['additional_processes']}
        if supplied != set(exclude_processes):
            raise ValueError(f'Explicit fixed-charge component must declare excluded source processes: {sorted(supplied)}')
        self.excluded_processes = sorted(supplied)
        self.lower = np.array([r['lower'] for r in rows], int)
        self.upper = np.array([r['upper'] for r in rows], int)
        self.temperature = self.data['temperatures_k']
        self.upsilon = np.array([r['upsilon'] for r in rows])
        if radiative_path is None:
            self.radiation = RadiativeLevels(dict(levels=self.data['levels'],
                lines=[r for r in rows if r['A_s1'] > 0], atomic_mass_u=39.948))
        else:
            self.radiation = RadiativeLevels.from_file(radiative_path)
            source_energy = np.array([r['energy_cm1'] for r in self.data['levels']])
            source_weight = np.array([r['degeneracy'] for r in self.data['levels']])
            if (len(source_energy) != len(self.radiation.energy) or
                    not np.array_equal(source_weight, self.radiation.weight) or
                    np.max(abs(source_energy-np.array([r['energy_cm1'] for r in self.radiation.record['levels']]))) > 1.):
                raise ValueError('Explicit radiation replacement requires matched ordered energies and weights')
        self.energy_ev = self.radiation.energy/E
        self.weight = self.radiation.weight
        self.delta_ev = self.energy_ev[self.upper]-self.energy_ev[self.lower]
        a0 = physical_constants['Bohr radius'][0]
        rydberg = physical_constants['Rydberg constant times hc in J'][0]
        self.rate_constant = np.sqrt(8/(np.pi*ME))*np.pi*a0*a0*rydberg/np.sqrt(KB)

    def rates(self, temperature_ev):
        """Excitation/de-excitation m^3/s, without extrapolating the source grid."""
        temperature = temperature_ev*E/KB
        if not np.isfinite(temperature) or not self.temperature[0] <= temperature <= self.temperature[-1]:
            raise ValueError('Electron temperature outside supplied collision-data range')
        index = min(np.searchsorted(self.temperature, temperature, side='right')-1, len(self.temperature)-2)
        fraction = np.log(temperature/self.temperature[index])/np.log(self.temperature[index+1]/self.temperature[index])
        upsilon = (1-fraction)*self.upsilon[:, index]+fraction*self.upsilon[:, index+1]
        common = self.rate_constant*upsilon/np.sqrt(temperature)
        return common/self.weight[self.lower]*np.exp(-self.delta_ev/temperature_ev), common/self.weight[self.upper]

    def generator(self, temperature_ev, electron_density_m3, *, radiation=True):
        if not np.isfinite(electron_density_m3) or electron_density_m3 < 0:
            raise ValueError('Finite nonnegative electron density required')
        up, down = [electron_density_m3*r for r in self.rates(temperature_ev)]
        count = len(self.energy_ev); matrix = np.zeros((count+2, count+2))
        for origin, destination, rates, energy_sign in ((self.lower, self.upper, up, -1), (self.upper, self.lower, down, 1)):
            np.add.at(matrix, (destination, origin), rates)
            np.add.at(matrix, (origin, origin), -rates)
            np.add.at(matrix[count], origin, energy_sign*self.delta_ev*rates)
        if radiation:
            r = self.radiation
            np.add.at(matrix, (r.lower, r.upper), r.A)
            np.add.at(matrix, (r.upper, r.upper), -r.A)
            np.add.at(matrix[count+1], r.upper, r.A*r.photon_energy/E)
        return matrix

    def step(self, state, dt, electron_density_m3, electrons_per_atom, *, radiation=True):
        """Second-order exponential midpoint with an independently evolved thermal reservoir."""
        state = np.asarray(state, float); count = len(self.energy_ev)
        if (state.shape != (count+2,) or not np.isfinite(state).all() or
                np.any(state[:count] < 0) or state[count] <= 0 or state[count+1] < 0 or
                not np.isfinite(dt) or dt <= 0 or not np.isfinite(electrons_per_atom) or electrons_per_atom <= 0):
            raise ValueError('Positive thermal energy, nonnegative populations/photons and positive step required')
        capacity = 1.5*electrons_per_atom
        temperature = state[count]/capacity
        half = expm(.5*dt*self.generator(temperature, electron_density_m3, radiation=radiation))@state
        middle_temperature = half[count]/capacity
        evolved = expm(dt*self.generator(middle_temperature, electron_density_m3, radiation=radiation))@state
        if not np.isfinite(evolved).all() or np.any(evolved[:count] < 0) or evolved[count] <= 0 or evolved[count+1] < 0:
            raise RuntimeError('Resolved excitation step left physical state domain; no clipping or energy repair')
        self.rates(evolved[count]/capacity)
        return evolved
