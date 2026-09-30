"""Bitwise source/state replay and independent higher-quadrature physical audit."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.face_radial_fluid import FaceRadialState
from .face_radial_study import initialize,STATE_KEYS,totals
from .face_radial_norms import FaceRadialSampler
from .face_radial_rate_study import calculate
from .radial_work_reference import RadialWorkReference
from .face_radial_records import provenance,sha


def audit_rates(directory,identity):
    results=[]
    for record in json.loads((directory/'rates.json').read_text()):
        path=directory/record['arrays_file'];assert sha(path)==record['arrays_sha256'] and record['source_sha256']==identity
        with np.load(path) as a:
            _,replay=calculate(record['spec'])
            assert set(a.files)==set(replay) and all(np.array_equal(a[k],v) for k,v in replay.items())
            ref,solver,_=initialize(record['spec'],RadialWorkReference())
            state=FaceRadialState(*(a[k] for k in STATE_KEYS));rate=FaceRadialState(*(a[k+'_rate'] for k in STATE_KEYS))
            metrics,energy=FaceRadialSampler(solver.mesh,record['spec'].get('degree',5),18).rates(state,rate,ref,0.)
            error=float(np.linalg.norm(energy-a['physical_energy_rate'])/np.linalg.norm(energy));assert error<5e-11
        results.append(dict(file=path.name,all_arrays_bitwise=True,independent_order18_norms=metrics,independent_energy_rate_relative=error))
        print(path.name,'rate audited',flush=True)
    return results


def audit_trajectories(directory,identity):
    results=[]
    for path in sorted(directory.glob('fluid_??.json')):
        record=json.loads(path.read_text());array_path=directory/record['arrays_file']
        assert sha(array_path)==record['arrays_sha256'] and record['source_sha256']==identity
        with np.load(array_path) as a:
            ref,solver,initial=initialize(record['spec'])
            def state(suffix,index=None):return FaceRadialState(*(a[k+suffix] if index is None else a[k+suffix][index] for k in STATE_KEYS))
            assert all(np.array_equal(getattr(initial,k),a[k+'_initial']) for k in STATE_KEYS)
            first_rate,_=solver.derivative(initial);assert all(np.array_equal(getattr(first_rate,k),a[k+'_initial_rate']) for k in STATE_KEYS)
            count=record['accepted_steps'];final=state('_final');caught=None
            if count:
                replay,_=solver.step(state('_before_final'),(count-1)*record['dt'],record['dt'])
                assert all(np.array_equal(getattr(replay,k),getattr(final,k)) for k in STATE_KEYS)
            if record['status']=='failed':
                before=[getattr(final,k).copy() for k in STATE_KEYS]
                try:solver.step(final,count*record['dt'],record['dt'])
                except (RuntimeError,ValueError) as error:caught=repr(error)
                assert caught==record['error'] and all(np.array_equal(v,getattr(final,k)) for v,k in zip(before,STATE_KEYS))
            grid=FaceRadialSampler(solver.mesh,record['spec'].get('degree',5),18);energy_errors=[];history_errors=[];norms=[]
            for i,t in enumerate(a['times']):
                s=state('_snapshots',i);metrics,energy=grid.measure(s,ref,float(t));norms.append(metrics)
                energy_errors.append(float(np.linalg.norm(energy-a['physical_energy_snapshots'][i])/np.linalg.norm(energy)))
                actual=np.array(totals(solver,s,float(t),a['spatial_ledger_snapshots'][i]));target=a['history'][int(round(t/record['dt']))]
                history_errors.append(float(np.max(abs(actual[:-2]-target[:-2])/np.maximum(abs(target[:-2]),1e-30))))
            final_norms,independent=grid.measure(final,ref,count*record['dt']);physical=solver.evaluate(final)[-1]
            final_error=float(np.linalg.norm(independent-physical)/np.linalg.norm(independent))
            assert max(history_errors)<1e-13
            if record['status']=='complete':assert max(energy_errors+[final_error])<2e-11
            h=a['history'];initial_energy=h[0,5]
            drift=(h[-1,5]-initial_energy)/initial_energy;spatial=h[-1,12]/initial_energy
            account=(h[-1,4]-h[0,4])/initial_energy
        results.append(dict(file=path.name,status=record['status'],initial_and_rate_bitwise=True,last_accepted_step_bitwise=bool(count),
                            rejected_step_replayed=caught is not None,max_independent_local_energy_relative=max(energy_errors+[final_error]),
                            quadrature_warning=max(energy_errors+[final_error])>=2e-11,max_history_relative=max(history_errors),
                            physical_energy_drift=drift,energy_account_drift=account,spatial_energy_defect=spatial,
                            finite_step_energy_defect=drift-account-spatial,independent_order18_norms=norms,final_norms=final_norms))
        print(path.name,'trajectory audited',record['status'],flush=True)
    return results


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args();d=args.directory
    identity=provenance()['source_sha256'];assert json.loads((d/'provenance.json').read_text())['source_sha256']==identity
    result=dict(source_sha256=identity,rates=audit_rates(d,identity) if (d/'rates.json').exists() else [],trajectories=audit_trajectories(d,identity))
    (d/'audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
