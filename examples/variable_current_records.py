"""Frozen source and raw-array identity for the variable-density component study."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import numpy as np
import scipy


def provenance():
    root=Path(__file__).resolve().parents[1];seen=set()
    def visit(module):
        path=root/(module.replace('.','/')+'.py')
        if module in seen or not path.is_file():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[a.name for a in node.names]
            else:continue
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for name in ('plasma_bench.__init__','examples.variable_current_verification',
                 'examples.variable_current_audit','tests.test_variable_compatible_current'):visit(name)
    files={m.replace('.','/')+'.py':hashlib.sha256((root/(m.replace('.','/')+'.py')).read_bytes()).hexdigest() for m in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
                python=sys.version,numpy=np.__version__,scipy=scipy.__version__)


def freeze(directory,record):
    root=Path(__file__).resolve().parents[1]
    with tarfile.open(directory/'sources.tar.gz','w:gz') as archive:
        for name in record['source_files']:
            archive.add(root/name,arcname=name)
            companion=(root/name).with_suffix('.md')
            if companion.exists():archive.add(companion,arcname=str(Path(name).with_suffix('.md')))
    (directory/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
