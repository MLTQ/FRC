"""Read personal LXCat text inputs and integrate their finite energy domains."""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import numpy as np
from scipy.constants import elementary_charge as EV, electron_mass as ME
from numpy.polynomial.legendre import leggauss


@dataclass
class CrossSection:
    kind: str
    target: str
    parameters: np.ndarray
    comments: tuple
    energy_ev: np.ndarray
    sigma_m2: np.ndarray


def read_cross_sections(path):
    """Read the explicit two-column LXCat format; retain source identity."""
    raw = Path(path).read_bytes(); lines = raw.decode('utf-8').splitlines()
    kinds = {'ELASTIC', 'EFFECTIVE', 'EXCITATION', 'IONIZATION', 'ATTACHMENT'}
    sections = []; i = 0
    while i < len(lines):
        kind = lines[i].strip(); i += 1
        if kind not in kinds:
            continue
        target = lines[i].strip(); i += 1
        parameters = np.array([]) if kind == 'ATTACHMENT' else np.array([float(v) for v in lines[i].split()])
        if kind != 'ATTACHMENT': i += 1
        comments = []
        while i < len(lines) and not re.fullmatch('-{5,}', lines[i].strip()):
            comments.append(lines[i]); i += 1
        if i == len(lines) or 'COLUMNS: Energy (eV) | Cross section (m2)' not in comments:
            raise ValueError('Explicit eV and m2 columns with table delimiters required')
        i += 1; rows = []
        while i < len(lines) and not re.fullmatch('-{5,}', lines[i].strip()):
            if lines[i].strip(): rows.append([float(v) for v in lines[i].split()])
            i += 1
        if i == len(lines): raise ValueError('Unterminated cross-section table')
        table = np.asarray(rows); i += 1
        if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2 or
                not np.isfinite(table).all() or np.any(table < 0) or
                np.any(np.diff(table[:, 0]) <= 0) or not np.isfinite(parameters).all()):
            raise ValueError('Finite nonnegative cross sections and increasing energies required')
        sections.append(CrossSection(kind, target, parameters, tuple(comments), table[:, 0], table[:, 1]))
    if not sections: raise ValueError('No supported cross sections found')
    return dict(input_path=str(Path(path).resolve()), input_sha256=hashlib.sha256(raw).hexdigest(), sections=sections)


class MaxwellianIntegral:
    """Integral of linear table segments only; no extrapolated high-energy tail."""
    def __init__(self, section, order=16):
        self.section = section; self.order = order
        nodes, weights = leggauss(order)
        e = section.energy_ev; s = section.sigma_m2
        fraction = .5*(nodes+1)
        self.energy = e[:-1, None]+np.diff(e)[:, None]*fraction
        sigma = s[:-1, None]+np.diff(s)[:, None]*fraction
        self.weighted = .5*np.diff(e)[:, None]*weights*sigma*self.energy
        self.factor = np.sqrt(8*EV/(np.pi*ME))

    def __call__(self, temperature_ev):
        t = float(temperature_ev)
        if not np.isfinite(t) or t <= 0: raise ValueError('Positive finite electron temperature required')
        return float(self.factor/t**1.5*np.sum(self.weighted*np.exp(-self.energy/t)))

    def electron_flux_fraction_above_table(self, temperature_ev):
        t = float(temperature_ev)
        if not np.isfinite(t) or t <= 0: raise ValueError('Positive finite electron temperature required')
        x = self.section.energy_ev[-1]/t
        return float((1+x)*np.exp(-x))
