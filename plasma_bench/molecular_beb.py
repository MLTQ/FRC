"""Attributed ground-molecule BEB ionization, without product-branch assumptions."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge, electron_mass
from scipy.special import gammaincc


BOHR_RADIUS_M = 5.29177210544e-11
RYDBERG_EV = 13.605693122990


@dataclass(frozen=True)
class MolecularBEB:
    gas: str
    orbitals: np.ndarray
    input_sha256: str
    maximum_energy_ev: float = 5000.

    @classmethod
    def from_nist(cls, directory, gas):
        if gas not in ('H2', 'N2', 'O2'):
            raise ValueError('Sourced ground-molecule data are available only for H2, N2 and O2')
        root = Path(directory)
        path = root/f'{gas}_orbitals.txt'
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        expected = json.loads((root/'manifest.json').read_text())['files'][path.name]['sha256']
        if digest != expected:
            raise ValueError('Orbital input does not match its recorded source hash')
        rows = []
        for line in path.read_text().splitlines()[3:]:
            fields = line.split('\t')
            if len(fields) != 8 or fields[4:] != ['1', 'No', 'none', 'none']:
                raise ValueError('Only the sourced Q=1, unmodified, single-ionization orbital formula is supported')
            rows.append([float(x) for x in fields[1:4]])
        orbitals = np.asarray(rows)
        if orbitals.ndim != 2 or orbitals.shape[1] != 3 or not np.isfinite(orbitals).all() or np.any(orbitals <= 0):
            raise ValueError('Positive finite binding, kinetic and occupation constants required')
        orbitals.setflags(write=False)
        return cls(gas, orbitals, digest)

    def orbital_cross_sections(self, energy_ev):
        energy = np.asarray(energy_ev, float)
        if not np.isfinite(energy).all() or np.any((energy < 0) | (energy > self.maximum_energy_ev)):
            raise ValueError('Energy must lie within the explicit 0–5000 eV comparison range')
        shape = energy.shape
        energy = energy.ravel()
        result = np.zeros((energy.size, len(self.orbitals)))
        for j, (binding, kinetic, occupation) in enumerate(self.orbitals):
            active = energy > binding
            delta = (energy[active]-binding)/binding
            t = 1+delta
            log_t = np.log1p(delta)
            factor = .5*log_t*delta*(t+1)/(t*t)+delta/t-log_t/(t+1)
            scale = 4*np.pi*BOHR_RADIUS_M**2*occupation*(RYDBERG_EV/binding)**2
            result[active, j] = scale*factor/(t+kinetic/binding+1)
        return result.reshape(shape+(len(self.orbitals),))

    def cross_section(self, energy_ev):
        return self.orbital_cross_sections(energy_ev).sum(axis=-1)

    def maxwellian_integral(self, temperature_ev, order=16):
        if not np.isfinite(temperature_ev) or temperature_ev <= 0:
            raise ValueError('Positive finite electron temperature required')
        if not isinstance(order, int) or order < 2:
            raise ValueError('Integer quadrature order >=2 required')
        nodes, weights = np.polynomial.legendre.leggauss(order)
        rate = 0.
        # Partition by excess energy over each orbital threshold. Exponential
        # scaling keeps low-temperature contributions visible to quadrature.
        for j, (binding, _, _) in enumerate(self.orbitals):
            limit = (self.maximum_energy_ev-binding)/temperature_ev
            cuts = [0.]
            edge = 1.
            while edge < limit:
                cuts.append(edge)
                edge *= 2
            cuts.append(limit)
            cuts = np.asarray(cuts)
            x = (cuts[:-1, None]+cuts[1:, None])/2+np.diff(cuts)[:, None]*nodes/2
            energy = binding+temperature_ev*x
            sigma = self.orbital_cross_sections(energy)[..., j]
            integral = np.sum(np.diff(cuts)[:, None]*weights/2*sigma*energy*np.exp(-x))
            rate += np.sqrt(8*elementary_charge/(np.pi*electron_mass*temperature_ev))*np.exp(-binding/temperature_ev)*integral
        return dict(rate_within_energy_range_m3_s=float(rate),
            maximum_energy_ev=self.maximum_energy_ev,
            electron_probability_above_range=float(gammaincc(1.5, self.maximum_energy_ev/temperature_ev)),
            temperature_ev=float(temperature_ev), quadrature_order=order)
