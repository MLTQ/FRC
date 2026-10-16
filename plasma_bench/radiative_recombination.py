"""Ground-state recombination sources with explicit escaping photon energy."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as EV


class RadiativeRecombination:
    """Total Maxwellian RR fits; capture-energy moment follows their derivative."""

    def __init__(self, record):
        self.record = record
        self.nuclear_charge = int(record['nuclear_charge'])
        self.rows = record['transitions']
        if [r['initial_charge'] for r in self.rows] != list(range(1, self.nuclear_charge+1)):
            raise ValueError('A complete ordered recombination ladder is required')
        for row in self.rows:
            if row['kind'] not in ('power_law', 'badnell'):
                raise ValueError('Unknown radiative recombination fit')
            if not np.isfinite(row['A_cm3_s']) or row['A_cm3_s'] <= 0:
                raise ValueError('Positive finite recombination prefactor required')
            fields = ('reference_temperature_ev',) if row['kind']=='power_law' else ('T0_ev','T1_ev')
            if any(not np.isfinite(row[k]) or row[k] <= 0 for k in fields):
                raise ValueError('Positive finite fit temperature scales required')
        self.range_known = np.array([r['minimum_temperature_ev'] is not None and
                                    r['maximum_temperature_ev'] is not None for r in self.rows])
        self.input_sha256 = None

    @classmethod
    def from_file(cls, path):
        raw = Path(path).read_bytes()
        model = cls(json.loads(raw))
        model.input_sha256 = hashlib.sha256(raw).hexdigest()
        return model

    def coefficients(self, temperature_j):
        t = float(temperature_j)/EV
        if not np.isfinite(t) or t <= 0:
            raise ValueError('Positive finite electron temperature required')
        rates = []; capture = []; within = []
        for row in self.rows:
            if row['kind']=='power_law':
                slope = -row['exponent']
                rate = row['A_cm3_s']*1e-6*(t/row['reference_temperature_ev'])**slope
            else:
                s0 = np.sqrt(t/row['T0_ev']); s1 = np.sqrt(t/row['T1_ev'])
                correction = row['C']*np.exp(-row['T2_ev']/t)
                b = row['B']+correction
                db = correction*row['T2_ev']/t
                rate = row['A_cm3_s']*1e-6/(s0*(1+s0)**(1-b)*(1+s1)**(1+b))
                slope = (-.5-.5*(1-b)*s0/(1+s0)-.5*(1+b)*s1/(1+s1)
                         +db*(np.log1p(s0)-np.log1p(s1)))
            mean_energy = temperature_j*(1.5+slope)
            if not np.isfinite([rate,mean_energy]).all() or rate < 0 or mean_energy < 0:
                raise ValueError('Rate fit has an inadmissible capture-energy moment; no clipping')
            rates.append(rate); capture.append(mean_energy)
            lo,hi = row['minimum_temperature_ev'],row['maximum_temperature_ev']
            within.append(lo is not None and hi is not None and lo <= t <= hi)
        return np.array(rates),np.array(capture),np.array(within)


class RecombiningPlasma:
    """Fixed-volume atomic source, assuming immediate cascade and photon escape."""

    def __init__(self, ionization, recombination, *, ionization_enabled=True):
        if ionization.nuclear_charge != recombination.nuclear_charge:
            raise ValueError('Ionization and recombination ladders must match')
        self.ionization = ionization
        self.recombination = recombination
        self.ionization_enabled = bool(ionization_enabled)
        self.charge = ionization.charge
        self.atomic_energy = ionization.atomic_energy

    def source(self, fractions, electron_energy_j_per_nucleus, nuclei_density_m3,
               *, heating_j_per_nucleus_s=0.):
        # Reuse the existing state-domain checks, without altering that component.
        df,de,di = self.ionization.source(fractions,electron_energy_j_per_nucleus,
            nuclei_density_m3,heating_j_per_nucleus_s=heating_j_per_nucleus_s)
        if not self.ionization_enabled:
            df = np.zeros_like(df); de = heating_j_per_nucleus_s
        f = np.asarray(fractions,float)
        temperature = di['temperature_j']
        if temperature is None:
            events = np.zeros(len(f)-1); capture = events.copy(); inside = np.zeros_like(events,bool)
        else:
            rates,capture,inside = self.recombination.coefficients(temperature)
            events = di['electron_density_m3']*f[1:]*rates
        df = df + np.r_[events,0.] - np.r_[0.,events]
        thermal_loss = float(events@capture)
        binding_release = float(events@self.ionization.ionization_energy)
        photons = thermal_loss+binding_release
        total = float(events.sum()); known = self.recombination.range_known
        diagnostic = dict(temperature_j=temperature,electron_density_m3=di['electron_density_m3'],
            recombination_events_per_nucleus_s=events,
            ionization_events_per_nucleus_s=di['events_per_nucleus_s'] if self.ionization_enabled else np.zeros_like(events),
            recombination_thermal_loss_j_per_nucleus_s=thermal_loss,
            recombination_binding_release_j_per_nucleus_s=binding_release,
            photon_power_j_per_nucleus_s=photons,
            unknown_rr_range_event_fraction=float(events[~known].sum()/total) if total else None,
            outside_known_rr_range_event_fraction=float(events[known & ~inside].sum()/total) if total else None,
            outside_ionization_fit_event_fraction=di['outside_fit_event_fraction'] if self.ionization_enabled else None)
        return df,de-thermal_loss,photons,diagnostic
