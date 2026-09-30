"""Bitwise and independent physical audits of all instantaneous coupling cases."""
import argparse
import json
from pathlib import Path
import numpy as np
from .regular_coupling_records import provenance,sha
from .regular_coupling_study import run_case,protocol
from .regular_coupling_reference import initialize
from .regular_coupling_norms import measure


def audit(directory):
    directory=Path(directory);identity=json.loads((directory/'provenance.json').read_text())
    if provenance()['source_sha256']!=identity['source_sha256']:raise ValueError('Current source mismatch')
    results=[]
    for i,spec in enumerate(protocol()):
        path=directory/f'case_{i:02d}.json';r=json.loads(path.read_text())
        assert r['source_sha256']==identity['source_sha256'] and r['status']=='complete' and r['spec']==spec
        raw=directory/r['arrays_file'];assert sha(raw)==r['arrays_sha256'];a=dict(np.load(raw))
        rr,aa=run_case(spec)
        for k,v in aa.items():assert np.array_equal(v,a[k]),(path.name,k)
        for k,v in rr.items():
            if k!='elapsed_seconds':assert v==r[k],(path.name,k)
        ref,d,u,m,_=initialize(spec);low=measure(d,m.mesh,ref,a,order=16);high=measure(d,m.mesh,ref,a,order=20)
        quadrature={k:abs(high[k]-low[k])/max(abs(high[k]),abs(low[k]),np.finfo(float).tiny) for k in low}
        current=float(np.sum(m.mesh.volumes*m.current_energy(a['current'])))
        grad=d.fields(u)['energy_gradient'];bulk=float(np.sum(d.volume[...,None]*grad*a['lorentz_rate']))
        kr=m.current_energy_rate(a['current'],a['current_rate'],a['fluid_rate'][...,0]);cr=m.current_energy_rate(a['current'],a['coupled_current_rate'],a['fluid_rate'][...,0])
        field_b=m.geometry.magnetic_samples(*m.geometry.unpack(a['magnetic']));field_b[...,2]+=ref.axial_bias
        def magnetic_rate(value):
            from scipy.constants import mu_0
            bt=m.geometry.magnetic_samples(*m.geometry.unpack(value))
            return float(np.sum(m.geometry.weight*np.sum(field_b*bt,axis=-1)/mu_0))
        prod=dict(current_energy=current,density_exchange=r['unresolved_density_exchange'],lorentz_work=bulk,
            current_energy_rate=float(np.sum(m.mesh.volumes*kr['total'])),magnetic_energy_rate=magnetic_rate(a['magnetic_rate']),
            coupled_current_energy_rate=float(np.sum(m.mesh.volumes*cr['total'])),coupled_magnetic_energy_rate=magnetic_rate(a['coupled_magnetic_rate']))
        diff={k:abs(high[k]-prod[k])/max(abs(high[k]),abs(prod[k]),np.finfo(float).tiny) for k in prod}
        result=dict(file=path.name,status='complete',arrays_and_metrics_bitwise=True,independent=high,quadrature_relative=quadrature,production_relative_difference=diff)
        results.append(result);print(path.name,'quadrature',max(quadrature.values()),'physical',max(diff.values()),flush=True)
    (directory/'audit.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],records=results),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');audit(p.parse_args().directory)
