"""Direct ionization charge source with independent atomic-energy accounting."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge


class ImpactIonization:
    def __init__(self, record):
        self.record = record
        rows = record['transitions']
        self.nuclear_charge = record['nuclear_charge']
        if [row['charge'] for row in rows] != list(range(self.nuclear_charge)):
            raise ValueError('A complete ordered ionization ladder is required')
        def column(name):
            return np.array([row[name] for row in rows], float)
        self.fit_energy = column('threshold_fit_ev') * elementary_charge
        self.ionization_energy = column('ionization_ev') * elementary_charge
        self.atomic_energy = np.r_[0., np.cumsum(self.ionization_energy)]
        self.prefactor = column('A_cm3_s') * 1e-6
        self.P = column('P'); self.X = column('X'); self.K = column('K')
        self.minimum_temperature = column('minimum_temperature_ev') * elementary_charge
        self.maximum_temperature = column('maximum_temperature_ev') * elementary_charge
        values = np.r_[self.fit_energy, self.ionization_energy, self.prefactor,
                       self.X, self.minimum_temperature, self.maximum_temperature]
        if not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError('Finite positive energies, prefactors and bounds required')
        if np.any(self.maximum_temperature < self.minimum_temperature):
            raise ValueError('Reversed rate-fit temperature interval')
        self.charge = np.arange(self.nuclear_charge + 1)
        self.input_sha256 = None

    @classmethod
    def from_file(cls, path):
        raw = Path(path).read_bytes()
        model = cls(json.loads(raw))
        model.input_sha256 = hashlib.sha256(raw).hexdigest()
        return model

    def coefficients(self, temperature_j):
        t = float(temperature_j)
        if not np.isfinite(t) or t <= 0:
            raise ValueError('Positive finite electron temperature required')
        u = self.fit_energy / t
        log_rate = (np.log(self.prefactor) + np.log1p(self.P*np.sqrt(u))
                    - np.log(self.X+u) + self.K*np.log(u) - u)
        return np.exp(log_rate), (t >= self.minimum_temperature) & (t <= self.maximum_temperature)

    def source(self, fractions, electron_energy_j_per_nucleus, nuclei_density_m3,
               *, heating_j_per_nucleus_s=0.):
        f = np.asarray(fractions, float)
        energy = float(electron_energy_j_per_nucleus)
        density = float(nuclei_density_m3)
        heating = float(heating_j_per_nucleus_s)
        if f.shape != self.charge.shape or not np.isfinite(f).all() or np.any(f < 0):
            raise ValueError('Finite nonnegative charge populations required; no clipping')
        if not np.isfinite([energy, density, heating]).all() or energy < 0 or density <= 0:
            raise ValueError('Finite energy/heating and positive nucleus density required')
        zbar = float(self.charge @ f)
        if zbar == 0:
            if energy != 0 or heating != 0:
                raise ValueError('Electron energy/heating requires a physical seed population')
            return np.zeros_like(f), 0., dict(electron_density_m3=0., temperature_j=None,
                events_per_nucleus_s=np.zeros(self.nuclear_charge), outside_fit_event_fraction=None,
                discrepant_threshold_event_fraction=None, ionization_power_j_per_nucleus_s=0.)
        temperature = energy / (1.5*zbar)
        rates, valid = self.coefficients(temperature)
        events = density*zbar*f[:-1]*rates
        df = np.r_[0., events] - np.r_[events, 0.]
        power = float(self.ionization_energy @ events)
        total = float(events.sum())
        discrepancy = abs(self.fit_energy/self.ionization_energy-1) > .05
        diagnostic = dict(electron_density_m3=density*zbar, temperature_j=temperature,
            events_per_nucleus_s=events, ionization_power_j_per_nucleus_s=power,
            outside_fit_event_fraction=float(events[~valid].sum()/total) if total else None,
            discrepant_threshold_event_fraction=float(events[discrepancy].sum()/total) if total else None)
        return df, heating-power, diagnostic
