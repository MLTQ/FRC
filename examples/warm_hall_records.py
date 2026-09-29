"""Frozen transitive source provenance and lossless warm Hall records."""
from pathlib import Path
import ast
import hashlib
import json
import sys
import numpy as np
import scipy
from .verification_records import write_record


def provenance():
    root = Path(__file__).resolve().parents[1]
    seen = set()
    def visit(module):
        path = root/(module.replace('.', '/')+'.py')
        if module in seen or not path.exists():
            return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and node.module:
                name = '.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module
                if name.startswith(('plasma_bench.', 'examples.')):
                    visit(name)
    for name in ('plasma_bench.__init__', 'examples.warm_hall_measurements', 'examples.warm_hall_records'):
        visit(name)
    files = {m.replace('.', '/')+'.py': hashlib.sha256((root/(m.replace('.', '/')+'.py')).read_bytes()).hexdigest() for m in sorted(seen)}
    return dict(source_files=files, source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__)


def save(directory, case, source, record, arrays):
    if provenance() != source:
        raise RuntimeError('Warm Hall measurement sources changed during run')
    path = directory/(case['name']+'.npz')
    np.savez_compressed(path, **arrays)
    result = dict(completed=True, configuration=case, provenance=source, **record,
                  fields=dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), units='SI; compatible magnetic energy, dual electron coefficients'))
    write_record(directory/(case['name']+'.json'), result)
    return result
