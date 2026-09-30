"""Source, independent physical and final-step audit of bulk thermal records."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.bulk_entropy_transport import BulkEntropyState
from .bulk_thermal_verification import calculate
from .bulk_thermal_evolution import initialize,totals
from .bulk_compression_reference import BulkCompressionReference
from .variable_current_reference import VariableCurrentReference
from .bulk_thermal_norms import rate_norms,transport_norms
from .bulk_thermal_records import provenance,sha


def common(directory,filename):
    live=provenance();frozen=json.loads((directory/'provenance.json').read_text())
    if live['source_sha256']!=frozen['source_sha256']:raise RuntimeError('Live source differs from frozen calculation')
    import hashlib
    with tarfile.open(directory/'sources.tar.gz') as archive:
        for name,wanted in frozen['source_files'].items():
            if hashlib.sha256(archive.extractfile(name).read()).hexdigest()!=wanted:raise RuntimeError('Source archive mismatch')
    records=json.loads((directory/filename).read_text())
    for r in records:
        if r['status']!='complete' or r['source_sha256']!=live['source_sha256'] or sha(directory/r['arrays_file'])!=r['arrays_sha256']:
            raise RuntimeError('Incomplete or modified raw record')
    return live,records


def audit_rates(directory):
    live,records=common(directory,'rates.json');cases=[]
    if len(records)!=14:raise RuntimeError('Expected fourteen rate cases')
    for index,r in enumerate(records):
        with np.load(directory/r['arrays_file']) as raw:arrays=dict(raw)
        _,replay=calculate(r['spec'])
        for key in arrays:np.testing.assert_array_equal(arrays[key],replay[key],err_msg=key)
        ref=VariableCurrentReference(**r['reference']);mesh=MHDMesh(ref.radius,ref.half_length,r['spec']['nr'],2*r['spec']['nr'])
        independent,local=rate_norms(mesh,ref,r['spec'],arrays,order=18)
        error={label:float(np.linalg.norm(local[label]-arrays[label+'_energy_rate'])/np.linalg.norm(local[label])) for label in ('electron','ion')}
        if max(error.values())>1e-7:raise RuntimeError('Independent local thermal derivative mismatch')
        cases.append(dict(spec=r['spec'],bitwise_rate_replay=True,independent_local_energy=error,norms=independent))
        print('audit rate',index,'complete',flush=True)
    result=dict(source_sha256=live['source_sha256'],records_sha256=sha(directory/'rates.json'),cases=cases)
    (directory/'audit.json').write_text(json.dumps(result,indent=2)+'\n')


def audit_transports(directory):
    live,records=common(directory,'transports.json');cases=[]
    if len(records)!=10:raise RuntimeError('Expected ten transport cases')
    for index,r in enumerate(records):
        with np.load(directory/r['arrays_file']) as raw:arrays=dict(raw)
        ref,solver,initial=initialize(r['spec'])
        for key in ('density','electron','ion'):np.testing.assert_array_equal(getattr(initial,key),arrays[key+'_initial'])
        before=BulkEntropyState(*(arrays[k+'_before_final'] for k in ('density','electron','ion')))
        last=solver.step(before,(r['steps']-1)*r['dt'],r['dt'])
        for key in ('density','electron','ion'):np.testing.assert_array_equal(getattr(last.state,key),arrays[key+'_final'])
        history=r['history']
        if len(history)!=r['steps']+1 or arrays['work_history'].shape!=(r['steps']+1,2):raise RuntimeError('Incomplete accepted history')
        for row in history:
            if not np.isfinite(row['particles']) or not row['min_density']>0:raise RuntimeError('Invalid accepted density')
            for label in ('electron','ion'):
                if not np.isfinite(list(row[label].values())).all() or row[label]['min_pressure']<=0:raise RuntimeError('Invalid accepted species pressure')
        np.testing.assert_allclose(arrays['work_history'][-1]-arrays['work_history'][-2],last.thermal_work,rtol=1e-10,atol=1e-22)
        norms=transport_norms(solver.mesh,ref,r['spec'],arrays,order=18);agreement=[]
        for sample,time in enumerate(arrays['times']):
            slot=int(round(time/r['dt']));row=norms[sample];saved=history[slot]
            diffs=dict(particles=abs(row['particles']/saved['particles']-1))
            for label in ('electron','ion'):
                for key in ('entropy','internal'):diffs[label+'_'+key]=abs(row[label][key]/saved[label][key]-1)
            if max(diffs.values())>1e-7:raise RuntimeError('Independent snapshot integral mismatch')
            agreement.append(diffs)
        cases.append(dict(spec=r['spec'],initial_and_last_step_bitwise=True,independent_integral_agreement=agreement,norms=norms))
        print('audit transport',index,'complete',flush=True)
    result=dict(source_sha256=live['source_sha256'],records_sha256=sha(directory/'transports.json'),cases=cases)
    (directory/'audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('rates',type=Path);parser.add_argument('transports',type=Path);args=parser.parse_args()
    audit_rates(args.rates);audit_transports(args.transports)
