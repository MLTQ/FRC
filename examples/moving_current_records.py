"""Frozen transitive source provenance for the moving-current reference study."""
from pathlib import Path
import ast
import hashlib
import json
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
            names=[]
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[n.name for n in node.names]
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for name in ('plasma_bench.__init__','plasma_bench.finite_mass_moments','examples.moving_current_measurements','examples.moving_current_records'):visit(name)
    files={m.replace('.','/')+'.py':hashlib.sha256((root/(m.replace('.','/')+'.py')).read_bytes()).hexdigest() for m in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__)
