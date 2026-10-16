"""Replace one declared ground-state ionization coefficient with a private table."""
import hashlib
from scipy.constants import elementary_charge as EV
from .lxcat_cross_sections import read_cross_sections, MaxwellianIntegral


class ArgonGroundCrossSection:
    """AtomicNetwork coefficient provider; all other charge stages use the base."""
    def __init__(self, base, path, *, order=16):
        self.base = base; self.data = read_cross_sections(path)
        rows = [r for r in self.data['sections'] if r.kind == 'IONIZATION' and r.target == 'Ar <-> Ar^+']
        if base.nuclear_charge != 18 or len(rows) != 1 or rows[0].parameters[0] <= 0:
            raise ValueError('Explicit unique ground-Ar ionization record and argon ladder required')
        self.rate = MaxwellianIntegral(rows[0], order)
        self.nuclear_charge = base.nuclear_charge
        self.charge = base.charge.copy(); self.atomic_energy = base.atomic_energy.copy()
        self.ionization_energy = base.ionization_energy.copy()
        self.fit_energy = base.fit_energy.copy(); self.fit_energy[0] = rows[0].parameters[0]*EV
        self.input_sha256 = hashlib.sha256((str(base.input_sha256)+self.data['input_sha256']).encode()).hexdigest()

    def coefficients(self, temperature_j):
        rates, inside = self.base.coefficients(temperature_j)
        rates[0] = self.rate(float(temperature_j)/EV)
        # This finite-energy table supplies no validated Maxwellian T interval.
        inside[0] = False
        return rates, inside
