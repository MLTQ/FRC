"""Transitive, separate provenance for the conservative finite-mass experiment."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy


def provenance():
    root=Path(__file__).resolve().parents[1];seen=set()
    def visit(module):
        path=root/(module.replace('.','/')+'.py')
        if module in seen or not path.exists():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[a.name for a in node.names]
            else:continue
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for module in ('plasma_bench.finite_mass_mhd','examples.finite_mass_rate_measurements','examples.finite_mass_evolution','examples.finite_mass_records'):visit(module)
    files={name.replace('.','/')+'.py':hashlib.sha256((root/(name.replace('.','/')+'.py')).read_bytes()).hexdigest() for name in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__)
