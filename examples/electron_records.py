"""Immutable source and lossless archive contracts for electron-pressure studies."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import scipy
from .verification_records import write_record


def provenance():
    paths = [f'plasma_bench/{s}.py' for s in ('dual_current', 'thermal_means', 'electron_pressure', 'electron_midpoint', 'mhd_mesh')]
    paths += [f'examples/{s}.py' for s in ('electron_reference', 'electron_cases', 'electron_records', 'electron_measurements', 'verification_records')]
    hashes = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
    return dict(source_files=hashes, source_sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__)


def relative(error, reference, weights):
    denominator = np.sum(weights*reference**2)
    return float(np.sqrt(np.sum(weights*error**2)/denominator)) if denominator else None


def save(directory, case, source, record, arrays):
    if provenance() != source:
        raise RuntimeError('Electron measurement source changed during run')
    path = directory/(case['name']+'.npz')
    np.savez_compressed(path, **arrays)
    result = dict(completed=True, configuration=case, provenance=source, **record,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; node/lumped electron state, cell g; not full fluid evolution'))
    write_record(directory/(case['name']+'.json'), result)
    return result
