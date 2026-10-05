"""Immutable complete source closure for bounded full moving verification."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy
from .radial_work_records import sha, freeze


def provenance():
    root = Path(__file__).resolve().parents[1]
    seen = set()
    def visit(module):
        path = root/(module.replace('.', '/')+'.py')
        if module in seen or not path.is_file():
            return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and node.module:
                names = ['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node, ast.Import):
                names = [v.name for v in node.names]
            else:
                continue
            for name in names:
                if name.startswith(('plasma_bench.', 'examples.', 'tests.')):
                    visit(name)
    for seed in ('plasma_bench.__init__', 'examples.__init__',
                 'examples.regular_plasma_protocol', 'tests.test_regular_plasma'):
        visit(seed)
    files = {m.replace('.', '/')+'.py': sha(root/(m.replace('.', '/')+'.py')) for m in sorted(seen)}
    return dict(source_files=files, source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__)
