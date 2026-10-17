"""Explicit-seed source closure for new implicit plasma calculations."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy
from .radial_work_records import sha, freeze


def provenance(seeds):
    root = Path(__file__).resolve().parents[1]; seen = set()
    def visit(module):
        path = root/(module.replace('.', '/')+'.py')
        if module in seen or not path.is_file(): return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and node.module:
                names = ['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node, ast.Import): names = [v.name for v in node.names]
            else: continue
            for name in names:
                if name.startswith(('plasma_bench.', 'examples.', 'tests.')): visit(name)
    for seed in ('plasma_bench.__init__', 'examples.__init__', __name__, *seeds): visit(seed)
    files = {m.replace('.', '/')+'.py': sha(root/(m.replace('.', '/')+'.py')) for m in sorted(seen)}
    runtime = {}
    if 'plasma_bench.precise_geometry' in seen:
        geometry = sys.modules.get('plasma_bench.precise_geometry')
        if geometry is not None:
            runtime.update(geometry_arithmetic=geometry.ARITHMETIC,
                           geometry_decimal_digits=geometry.DIGITS,
                           geometry_precision_bits=geometry.BITS)
            if geometry.ARITHMETIC == 'mpfr':
                runtime.update(gmpy2=geometry.gmpy2.version(),
                               mpfr=geometry.gmpy2.mpfr_version(),
                               gmp=geometry.gmpy2.mp_version())
    return dict(source_files=files, source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__, **runtime)
