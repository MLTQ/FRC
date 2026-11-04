"""Complete electrostatic restart state, with physical/source compatibility checks."""
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
import numba

from .iec_particles import Particles


PHYSICS = ('iec_field', 'iec_particles', 'iec_mesh', 'iec_vacuum', 'hops_drive', 'iec_diagnostics')
SCALARS = ('initial', 'loss_energy', 'injected_energy', 'source_work', 'injected_charge',
           'time_s', 'energy', 'space', 'cage_voltage')
ARRAYS = ('phi', 'charge', 'electric', 'electrodes', 'lost_charge', 'lost_count')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def physics_hashes():
    root = Path(__file__).parent
    return {name:digest(root/(name+'.py')) for name in PHYSICS}


def libraries():
    return dict(numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,
                machine=platform.machine())


def write_checkpoint(path, config, step, dt, runtime):
    """Atomic full-step checkpoint, after the second kick and saved diagnostics."""
    p, drive, rng = (runtime[k] for k in ('p', 'drive', 'rng'))
    meta = dict(version=1, config=config, step=step, dt=dt,libraries=libraries(),
                weight=float(p.weight), ion_mass=float(p.ion_mass),
                scalars={k:float(runtime[k]) for k in SCALARS},
                rng=rng.bit_generator.state,
                drive={k:v for k,v in vars(drive).items() if k != 'config'},
                history_rows=len(runtime['history']), history=runtime['history'],
                elapsed_s=runtime['history'][-1]['elapsed_s'],
                mesh_sha256=runtime['mesh_sha256'], physics_hashes=runtime['physics_hashes'])
    values={k:runtime[k] for k in ARRAYS}
    values.update({k:getattr(p,k) for k in ('x','v','cells','species')})
    values['metadata_json'] = np.array(json.dumps(meta,allow_nan=False))
    path=Path(path); temporary=path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**values)
    temporary.replace(path)


def load_checkpoint(path, config, mesh_path):
    with np.load(path,allow_pickle=False) as archive:
        meta=json.loads(str(archive['metadata_json']))
        values={k:archive[k].copy() for k in archive.files if k != 'metadata_json'}
    if meta['version'] != 1 or meta['config'] != config:
        raise ValueError('Checkpoint configuration mismatch')
    if meta['libraries'] != libraries():
        raise ValueError('Checkpoint numerical libraries or architecture differ')
    if meta['mesh_sha256'] != digest(mesh_path) or meta['physics_hashes'] != physics_hashes():
        raise ValueError('Checkpoint geometry or physics source mismatch')
    if not all(np.isfinite(v).all() for v in values.values()):
        raise ValueError('Nonfinite checkpoint')
    p=Particles(*(values.pop(k) for k in ('x','v','cells','species')),meta['weight'],meta['ion_mass'])
    return meta,values,p
