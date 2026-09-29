"""Shared source/configuration and lossless archive contracts for Hall studies."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import scipy
from .verification_records import write_record


def provenance():
    production = ('compatible_magnetic', 'compatible_poynting', 'hall_curl', 'hall_induction',
                  'hall_compatible', 'hall_matrices', 'hall_midpoint', 'hall_fluid', 'magnetic_quadrature',
                  'toroidal_state', 'toroidal_mhd', 'toroidal_flux', 'toroidal_reconstruction', 'toroidal_induction',
                  'mhd', 'mhd_mesh', 'mhd_state', 'mhd_flux', 'mhd_induction', 'mhd_reconstruction', 'mhd_positivity')
    examples = ('compatible_hall_references', 'compatible_hall_cases', 'compatible_hall_records',
                'compatible_hall_spatial', 'compatible_hall_evolution', 'compatible_hall_measurements',
                'compressible_hall_reference', 'hall_manufactured', 'hall_work_snapshots', 'verification_records')
    names = [f'plasma_bench/{name}.py' for name in production]+[f'examples/{name}.py' for name in examples]
    files = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}
    return dict(source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(), source_files=files,
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__)


def relative(error, reference, weight):
    denominator = np.sum(weight*reference*reference)
    return float(np.sqrt(np.sum(weight*error*error)/denominator)) if denominator > 0 else None


def save(directory, case, record, arrays, source):
    path = directory/(case['name']+'.npz')
    np.savez_compressed(path, **arrays)
    if provenance() != source:
        raise RuntimeError('Compatible Hall sources changed during measurement')
    result = dict(completed=True, configuration=case, provenance=source, **record,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; explicit compatible magnetic energy; see module companions'))
    write_record(directory/(case['name']+'.json'), result)
    return result
