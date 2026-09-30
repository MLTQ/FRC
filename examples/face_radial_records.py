"""Source closure and atomic saved records for the common radial-face study."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy
from .radial_work_records import sha,freeze


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
    for seed in ('examples.face_radial_study','examples.face_radial_rate_study',
                 'examples.face_radial_audit','tests.test_face_radial_fluid'):visit(seed)
    files={m.replace('.','/')+'.py':sha(root/(m.replace('.','/')+'.py')) for m in sorted(seen)}
    return dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
                python=sys.version,numpy=np.__version__,scipy=scipy.__version__)


def save(directory,name,record,arrays):
    """Array first, then atomic metadata publication; previous checkpoint survives."""
    suffix=f"_step{record['accepted_steps']:06d}" if record['status']=='running' else ''
    path=directory/(name+suffix+'.npz');temporary=path.with_suffix('.tmp')
    with temporary.open('wb') as stream:np.savez_compressed(stream,**arrays)
    temporary.replace(path)
    record=dict(record,arrays_file=path.name,arrays_sha256=sha(path))
    target=directory/(name+('_checkpoint' if suffix else '')+'.json')
    temporary=target.with_suffix('.tmp');temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(target)
    return record
