"""Explicit local ADF07 ionization rates with their source state identifiers."""
import hashlib
from pathlib import Path
import re
import numpy as np


def read_ionization(path):
    raw = Path(path).read_bytes()
    lines = raw.decode('ascii').splitlines()
    count = int(lines[0].split()[0]); position = 1; rows = []
    pattern = (r'\s*([A-Za-z]+)\+\s*(\d+)/([A-Za-z]+)\+\s*(\d+)/\s*(\d+)/'
               r'I\.P\.\s*=\s*([\d.]+)/ICODE\s*=\s*(\d+)/FCODE\s*=\s*(\d+)/ISEL\s*=\s*(\d+)')
    for _ in range(count):
        match = re.fullmatch(pattern, lines[position].strip())
        if match is None:
            raise ValueError(f'Unsupported ADF07 reaction header on line {position+1}')
        element, charge, final_element, final_charge, nt, potential, initial, final, index = match.groups()
        position += 1; nt = int(nt); values = []
        while len(values) < 2*nt:
            values.extend(float(value.replace('D', 'E')) for value in lines[position].split())
            position += 1
        temperature, rate = np.asarray(values[:nt]), np.asarray(values[nt:])*1e-6
        if (len(values) != 2*nt or nt < 2 or element != final_element or
                int(final_charge) != int(charge)+1 or not np.isfinite(values).all() or
                np.any(temperature <= 0) or np.any(np.diff(temperature) <= 0) or
                np.any(rate <= 0) or float(potential) <= 0):
            raise ValueError('Positive ordered temperatures/rates and adjacent charge transition required')
        rows.append(dict(element=element, charge=int(charge), final_charge=int(final_charge),
            initial_metastable=int(initial), final_metastable=int(final), index=int(index),
            effective_potential_cm1=float(potential), temperature_ev=temperature, rate_m3_s=rate))
    if len({row['index'] for row in rows}) != count:
        raise ValueError('Unique source reaction indices required')
    return dict(rows=rows, input_path=str(Path(path).resolve()),
                input_sha256=hashlib.sha256(raw).hexdigest())


def coefficient(row, temperature_ev):
    """Positive log-log interpolation inside the supplied interval, no extrapolation."""
    grid = row['temperature_ev']
    if not np.isfinite(temperature_ev) or not grid[0] <= temperature_ev <= grid[-1]:
        raise ValueError('Temperature outside supplied ADF07 ionization interval')
    return float(np.exp(np.interp(np.log(temperature_ev), np.log(grid), np.log(row['rate_m3_s']))))
