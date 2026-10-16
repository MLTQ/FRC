"""Resolved radiative decays, stored level energy and line absorption."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import atomic_mass, c, h
from scipy.special import voigt_profile


class RadiativeLevels:
    """A finite atomic level set with explicit spontaneous photon production."""

    def __init__(self, record):
        self.record = record
        levels, lines = record['levels'], record['lines']
        self.energy = np.array([x['energy_cm1'] for x in levels])*100*h*c
        self.weight = np.array([x['degeneracy'] for x in levels], float)
        self.lower = np.array([x['lower'] for x in lines], int)
        self.upper = np.array([x['upper'] for x in lines], int)
        self.A = np.array([x['A_s1'] for x in lines], float)
        self.mass = record['atomic_mass_u']*atomic_mass
        if (not np.isfinite(self.energy).all() or self.energy[0] != 0
                or np.any(np.diff(self.energy) <= 0) or np.any(self.weight <= 0)
                or np.any(self.lower < 0) or np.any(self.upper >= len(levels))
                or np.any(self.upper <= self.lower) or not np.isfinite(self.A).all()
                or np.any(self.A <= 0) or self.mass <= 0):
            raise ValueError('Ordered finite levels and positive physical transitions required')
        self.photon_energy = self.energy[self.upper]-self.energy[self.lower]
        # Energy differences define photon energies and wavelengths consistently.
        self.frequency = self.photon_energy/h
        self.wavelength = c/self.frequency
        self.decay_rate = np.bincount(self.upper, weights=self.A, minlength=len(levels))
        self.input_sha256 = None

    @classmethod
    def from_file(cls, path):
        raw = Path(path).read_bytes()
        model = cls(json.loads(raw))
        model.input_sha256 = hashlib.sha256(raw).hexdigest()
        return model

    def source(self, populations):
        """Return level transfers and produced photons per line, per atom / s."""
        f = np.asarray(populations, float)
        if f.shape != self.energy.shape or not np.isfinite(f).all() or np.any(f < 0):
            raise ValueError('Finite nonnegative level populations required; no repair')
        events = self.A*f[self.upper]
        change = (np.bincount(self.lower, weights=events, minlength=len(f))
                  - np.bincount(self.upper, weights=events, minlength=len(f)))
        return change, events

    def widths(self, temperature_j):
        """Gaussian frequency standard deviation and natural Lorentz HWHM [Hz]."""
        if not np.isfinite(temperature_j) or temperature_j <= 0:
            raise ValueError('Positive atom temperature k_B T in joules required')
        sigma = self.frequency/c*np.sqrt(temperature_j/self.mass)
        gamma = (self.decay_rate[self.upper]+self.decay_rate[self.lower])/(4*np.pi)
        return sigma, gamma

    def absorption(self, level_density_m3, temperature_j, frequency_offset_hz=0.):
        """Net absorption [m^-1], including stimulated emission; negative is gain."""
        n = np.asarray(level_density_m3, float)
        if n.shape != self.energy.shape or not np.isfinite(n).all() or np.any(n < 0):
            raise ValueError('Finite nonnegative level densities required')
        sigma, gamma = self.widths(temperature_j)
        profile = voigt_profile(np.asarray(frequency_offset_hz), sigma, gamma)
        inversion = n[self.lower]*self.weight[self.upper]/self.weight[self.lower]-n[self.upper]
        return self.A*self.wavelength**2/(8*np.pi)*inversion*profile
