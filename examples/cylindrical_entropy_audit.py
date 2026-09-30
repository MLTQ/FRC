"""Hash, original-equation, physical-energy and accepted-state replay audit."""
import argparse
import json
from pathlib import Path
import numpy as np
from .cylindrical_entropy_records import provenance,sha
from .cylindrical_entropy_study import initialize,rate_case,history_row
from .cylindrical_entropy_spectrum import calculate
from .cylindrical_entropy_reference import physical,rates,measure


def audit(directory):
    identity=provenance();saved=json.loads((directory/'provenance.json').read_text())
    if saved['source_sha256']!=identity['source_sha256']:raise ValueError('Source identity mismatch')
    result=[]
    for path in sorted(directory.glob('*.json')):
        if path.stem in ('audit','provenance') or 'checkpoint' in path.stem:continue
        record=json.loads(path.read_text())
        if record.get('source_sha256')!=identity['source_sha256']:raise ValueError(f'Wrong record source {path}')
        raw=directory/record['arrays_file']
        if sha(raw)!=record['arrays_sha256']:raise ValueError(f'Raw-array hash mismatch {path}')
        arrays=dict(np.load(raw));entry=dict(file=path.name)
        if path.stem.startswith(('rates_','spectra_')):
            current,replay=(rate_case(record['spec']) if path.stem.startswith('rates_') else calculate(record['spec']))
            for key,value in replay.items():
                if not np.array_equal(value,arrays[key]):raise ValueError(f'Non-bitwise replay {path}:{key}')
            for key,value in current.items():
                if value!=record[key]:raise ValueError(f'Metric replay mismatch {path}:{key}')
            entry['replayed_arrays']=len(replay)
        else:
            ref,d,initial=initialize(record['spec']);u=arrays['state'];du,diagnostic=d.derivative(u)
            if not np.array_equal(initial,arrays['initial']) or not np.array_equal(du,arrays['rate']):raise ValueError(f'Initial/rate replay mismatch {path}')
            if record['status']=='complete':
                if record['accepted_steps']!=record['steps']:raise ValueError('Incomplete completed trajectory')
                if not np.array_equal(d.step(arrays['previous'],record['dt']),u):raise ValueError(f'Last step mismatch {path}')
            elif record['status']=='failed':
                before=u.copy()
                try:d.step(u,record['dt'])
                except ValueError:pass
                else:raise ValueError('Failure did not replay')
                if not np.array_equal(before,u):raise ValueError('Rejected step mutated state')
            else:raise ValueError('Nonterminal record passed as final')
            for state,step,row in zip(arrays['snapshots'],arrays['snapshot_steps'],arrays['history'],strict=True):
                if not np.array_equal(history_row(d,state,ref,step*record['dt']),row):raise ValueError(f'Physical history mismatch {path}')
            independent=rates(u,du,d.r,ref)['energy'];scale=max(np.linalg.norm(independent),1e-300)
            entry['independent_local_energy_rate_relative']=float(np.linalg.norm(independent-diagnostic['physical_energy_rate'])/scale)
            for key,order,factor,time,state in (('initial_metrics',32,4,0.,initial),('final_metrics',32,4,record['duration'],u),('higher_metrics',48,8,record['duration'],u)):
                if measure(d,state,ref,time,order,factor)!=record[key]:raise ValueError(f'Physical integration mismatch {path}:{key}')
            lo=record['final_metrics']['continuous_energy'];hi=record['higher_metrics']['continuous_energy']
            entry.update(status=record['status'],history_states=len(arrays['snapshots']),higher_quadrature_relative_difference=abs(lo-hi)/abs(hi))
        result.append(entry)
    output=dict(source_sha256=identity['source_sha256'],records=result)
    (directory/'audit.json').write_text(json.dumps(output,indent=2)+'\n')
    return output


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args()
    print(json.dumps(audit(args.directory),indent=2))


if __name__=='__main__':main()
