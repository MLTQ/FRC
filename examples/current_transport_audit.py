"""Bitwise replay and independent physical integration of transport evidence."""
import argparse
import json
from pathlib import Path
import numpy as np
from .current_transport_records import provenance,sha
from .current_transport_study import run_case,protocol
from .current_transport_norms import measure
from .regular_coupling_reference import initialize


def audit(directory):
    directory=Path(directory);identity=json.loads((directory/'provenance.json').read_text())
    if provenance()['source_sha256']!=identity['source_sha256']:raise ValueError('Current source mismatch')
    results=[]
    for i,spec in enumerate(protocol()):
        path=directory/f'case_{i:02d}.json';r=json.loads(path.read_text())
        assert r['source_sha256']==identity['source_sha256'] and r['status']=='complete' and r['spec']==spec
        raw=directory/r['arrays_file'];assert sha(raw)==r['arrays_sha256'];a=dict(np.load(raw));rr,aa=run_case(spec)
        assert set(a)==set(aa)
        for k,v in aa.items():assert np.array_equal(v,a[k]),(path.name,k)
        for k,v in rr.items():
            if k!='elapsed_seconds':assert v==r[k],(path.name,k)
        ref,d,u,m,field=initialize(spec);low=measure(d,m.mesh,ref,a,m.delta,order=16);high=measure(d,m.mesh,ref,a,m.delta,order=20)
        scale=high['powers']['continuous_rate_scale']
        quadrature={k:dict(absolute=abs(high['powers'][k]-v),relative=abs(high['powers'][k]-v)/max(abs(high['powers'][k]),abs(v),np.finfo(float).tiny)) for k,v in low['powers'].items()}
        production={k:r[k] for k in ('advection_work','stretching_work','self_drift_work','density_exchange')}
        kr=m.current_energy_rate(a['current'],a['current_rate'],a['fluid_rate'][...,0]);br=m.magnetic_energy_rate(field,a['magnetic_rate'])
        production.update(current_energy_rate=float(np.sum(m.mesh.volumes*kr['total'])),magnetic_energy_rate=float(np.sum(m.mesh.volumes*br)))
        differences={k:dict(absolute=abs(high['powers'][k]-v),relative=abs(high['powers'][k]-v)/max(abs(high['powers'][k]),abs(v),np.finfo(float).tiny),scaled_to_physical_rate=abs(high['powers'][k]-v)/scale) for k,v in production.items()}
        result=dict(file=path.name,arrays_and_metrics_bitwise=True,independent=high,quadrature_power=quadrature,production_power_difference=differences,
            grid_stress_work=r['stress_work'],continuous_stress_work_minus_grid=high['powers']['stress_rate_physical_work']-r['stress_work'],grid_total_rate=r['total_energy_rate'])
        results.append(result)
        (directory/'audit_progress.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],terminal=False,records=results),indent=2)+'\n')
        print(path.name,'independent force',high['norms']['stress_rate']['relative_L2'],'continuous energy',high['powers']['continuous_total_rate'],flush=True)
    assert provenance()['source_sha256']==identity['source_sha256']
    (directory/'audit.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],terminal=True,records=results),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');audit(p.parse_args().directory)
