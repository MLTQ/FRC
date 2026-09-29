"""Source identities and lossless records for the moving species study."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import scipy
from .verification_records import write_record


def provenance():
    production = ('dual_current', 'dual_mass', 'bulk_entropy', 'two_temperature', 'warm_mhd',
                  'toroidal_mhd', 'toroidal_state', 'toroidal_flux', 'toroidal_reconstruction', 'toroidal_induction',
                  'mhd', 'mhd_mesh', 'mhd_state', 'mhd_flux', 'mhd_induction', 'mhd_reconstruction', 'mhd_positivity',
                  'magnetic_quadrature', 'compatible_magnetic', 'magnetic_dissipation')
    examples = ('warm_reference', 'warm_cases', 'warm_records', 'warm_measurements', 'verification_records')
    paths = [f'plasma_bench/{s}.py' for s in production]+[f'examples/{s}.py' for s in examples]
    hashes = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
    return dict(source_files=hashes, source_sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__)


def save(directory, case, source, record, arrays):
    if provenance() != source:
        raise RuntimeError('Warm source changed during measurement')
    path = directory/(case['name']+'.npz')
    np.savez_compressed(path, **arrays)
    result = dict(completed=True, configuration=case, provenance=source, **record,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; primal conservative state and dual electron coefficients'))
    write_record(directory/(case['name']+'.json'), result)
    return result
