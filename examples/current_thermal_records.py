"""Frozen transitive current/thermal source identity and strict lossless records."""
from pathlib import Path
import ast
import hashlib
import json
import sys
import numpy as np
import scipy
from .verification_records import write_record


def provenance():
    root=Path(__file__).resolve().parents[1];seen=set()
    def visit(module):
        path=root/(module.replace('.','/')+'.py')
        if module in seen or not path.exists():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            names=[]
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[n.name for n in node.names]
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for name in ('plasma_bench.__init__','examples.current_thermal_measurements','examples.current_thermal_records'):visit(name)
    files={m.replace('.','/')+'.py':hashlib.sha256((root/(m.replace('.','/')+'.py')).read_bytes()).hexdigest() for m in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
                python=sys.version,numpy=np.__version__,scipy=scipy.__version__)


def save(directory,case,source,record,arrays):
    if provenance()!=source:raise RuntimeError('Current/thermal source changed during measurement')
    path=directory/(case['name']+'.npz');np.savez_compressed(path,**arrays)
    result=dict(configuration=case,provenance=source,**record,
        fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),units='SI; face powers W; entropy U=pressure/(n/nref)^(gamma-1)'))
    write_record(directory/(case['name']+'.json'),result)
    return result
