"""Read an explicitly supplied Maxwellian ADF04 excitation file, without downloading."""
import hashlib
from pathlib import Path
import re
import numpy as np


def adas_number(text):
    """ADAS omits E in numbers such as 5.00+02."""
    try:
        return float(text)
    except ValueError:
        return float(re.sub(r'(?<=\d)([+-]\d+)$', r'e\1', text))


def read_levels(path):
    raw = Path(path).read_bytes(); lines = raw.decode('ascii').splitlines()
    header = re.match(r'\s*([A-Za-z]+)\+\s*(\d+)\s+(\d+)\s+(\d+)\s+', lines[0])
    if header is None:
        raise ValueError('Expected a resolved-ion ADF04 header')
    levels = []; i = 1
    pattern = r'\s*(\d+)\s+(.*?)\s+\((\d+)\)(\d+)\(\s*([\d.]+)\)\s*([\d.]+)'
    while lines[i].strip() != '-1':
        match = re.match(pattern, lines[i])
        if match is None:
            raise ValueError(f'Unsupported level record on line {i+1}')
        index, name, multiplicity, orbital, xj, energy = match.groups()
        if int(index) != len(levels)+1:
            raise ValueError('Contiguous source level ordering required')
        levels.append(dict(index=int(index), configuration=name.strip(), multiplicity=int(multiplicity),
            orbital_L=int(orbital), source_xj=float(xj), degeneracy=2*float(xj)+1, energy_cm1=float(energy)))
        i += 1
    i += 1; values = lines[i].split()
    if int(values[1]) != 3:
        raise ValueError('Only type-3 Maxwellian effective collision strengths are supported')
    temperatures = np.array([adas_number(v) for v in values[2:]])
    if not np.isfinite(temperatures).all() or np.any(temperatures <= 0) or np.any(np.diff(temperatures) <= 0):
        raise ValueError('Increasing positive source temperatures required')
    i += 1; transitions = []; additional_processes = []
    while not lines[i].lstrip().startswith('-1'):
        if not lines[i].startswith(' '):
            additional_processes.append(dict(code=lines[i][0], source_line=i+1, record=lines[i]))
            i += 1
            continue
        parts = lines[i].split(); upper, lower = map(int, parts[:2])
        values = [adas_number(v) for v in parts[2:]]
        if len(values) not in (len(temperatures)+1, len(temperatures)+2):
            raise ValueError('Unexpected excitation/Born-limit columns')
        if not 1 <= lower < upper <= len(levels) or not np.isfinite(values).all() or min(values[:len(temperatures)+1]) < 0:
            raise ValueError('Ordered valid transitions and nonnegative literal source data required')
        transitions.append(dict(lower=lower-1, upper=upper-1, A_s1=values[0],
            upsilon=values[1:len(temperatures)+1], born_limit=values[-1] if len(values)==len(temperatures)+2 else None))
        i += 1
    pairs = [(r['lower'], r['upper']) for r in transitions]
    if len(set(pairs)) != len(pairs):
        raise ValueError('Duplicate excitation rows require explicit physical interpretation')
    return dict(element=header[1], charge=int(header[2]), atomic_number=int(header[3]),
        levels=levels, temperatures_k=temperatures, transitions=transitions,
        additional_processes=additional_processes,
        input_path=str(Path(path).resolve()), input_sha256=hashlib.sha256(raw).hexdigest(),
        interpretation='Source xj gives (g-1)/2; it need not be a true J for LS term-bundled data.')
