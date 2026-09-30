"""Bitwise rate replay and independent high-quadrature component diagnosis."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.log_fluid import LogFluidState
from .radial_work_rate_study import calculate
from .radial_work_study import initialize,STATE_KEYS
from .radial_work_reference import RadialWorkReference
from .log_fluid_norms import FluidSampler
from .radial_work_records import provenance,sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    records=json.loads((d/'rates.json').read_text());identity=provenance();audits=[]
    for record in records:
        assert record['source_sha256']==identity['source_sha256'];path=d/record['arrays_file']
        assert sha(path)==record['arrays_sha256'];a=np.load(path);_,replay=calculate(record['spec'])
        assert set(a.files)==set(replay) and all(np.array_equal(a[k],v) for k,v in replay.items())
        _,solver,_=initialize(dict(kind='contact',**record['spec']));ref=RadialWorkReference()
        state=LogFluidState(*(a[k] for k in STATE_KEYS));rate=LogFluidState(*(a[k+'_rate'] for k in STATE_KEYS))
        grid=FluidSampler(solver.mesh,record['spec'].get('degree',5),18);metrics,energy=grid.rates(state,rate,ref,0.)
        f=grid.physical(state,ref);nt=grid.density(rate.density);mt=grid.momentum(rate.momentum)
        vt=(mt-f['velocity']*nt[...,None])/f['n'][...,None];target=ref.rates(grid.r,grid.z)['velocity']
        components={key:grid.metric(vt[...,k],target[...,k]) for k,key in enumerate(('radial','azimuthal','axial'))}
        agreement=float(np.linalg.norm(energy-a['physical_energy_rate'])/np.linalg.norm(energy))
        assert agreement<5e-11
        audits.append(dict(file=path.name,all_arrays_bitwise=True,higher_quadrature=metrics,velocity_components=components,
                           independent_energy_rate_relative=agreement))
        print(path.name,'audited',components['radial']['axis_relative_L2'],flush=True)
    source=Path(__file__)
    (d/'rate_audit.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],audit_sha256=sha(source),records=audits),indent=2)+'\n')
    (d/'audit_source.py').write_bytes(source.read_bytes());(d/'audit_source.md').write_bytes(source.with_suffix('.md').read_bytes())


if __name__=='__main__':main()
