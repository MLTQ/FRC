"""Immutable source and raw evidence for conservative bulk thermal transport."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import numpy as np
import scipy


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance():
    root=Path(__file__).resolve().parents[1];seen=set()
    def visit(module):
        path=root/(module.replace('.','/')+'.py')
        if module in seen or not path.is_file():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[v.name for v in node.names]
            else:continue
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for seed in ('plasma_bench.__init__','examples.bulk_thermal_verification','examples.bulk_thermal_evolution',
                 'examples.bulk_thermal_audit','tests.test_compatible_bulk_thermal','tests.test_bulk_entropy_transport',
                 'tests.test_normal_momentum_faces'):visit(seed)
    files={m.replace('.','/')+'.py':sha(root/(m.replace('.','/')+'.py')) for m in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
                python=sys.version,numpy=np.__version__,scipy=scipy.__version__)


def freeze(directory,identity):
    root=Path(__file__).resolve().parents[1]
    with tarfile.open(directory/'sources.tar.gz','w:gz') as archive:
        for name in identity['source_files']:
            path=root/name;archive.add(path,arcname=name)
            companion=path.with_suffix('.md')
            if companion.exists():archive.add(companion,arcname=str(Path(name).with_suffix('.md')))
    (directory/'provenance.json').write_text(json.dumps(identity,indent=2)+'\n')
