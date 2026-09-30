"""Replay stored fluid steps and independently integrate all saved physical states."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.log_fluid import LogFluidState
from .log_fluid_study import initialize,STATE_KEYS,totals
from .log_fluid_norms import FluidSampler
from .log_fluid_records import provenance,sha


def audit(directory):
    identity=json.loads((directory/'provenance.json').read_text())
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Audit source mismatch')
    records=[]
    for path in sorted(directory.glob('fluid_*.json')):
        record=json.loads(path.read_text());arrays=directory/record['arrays_file']
        assert sha(arrays)==record['arrays_sha256'];assert record['source_sha256']==identity['source_sha256']
        a=np.load(arrays);ref,solver,initial=initialize(record['spec'])
        def state(suffix,index=None):return LogFluidState(*(a[k+suffix] if index is None else a[k+suffix][index] for k in STATE_KEYS))
        assert all(np.array_equal(getattr(initial,k),a[k+'_initial']) for k in STATE_KEYS)
        first_rate,_=solver.derivative(initial)
        assert all(np.array_equal(getattr(first_rate,k),a[k+'_initial_rate']) for k in STATE_KEYS)
        if record['status']=='complete':
            replay,d=solver.step(state('_before_final'),record['end_time']-record['dt'],record['dt'])
            assert all(np.array_equal(getattr(replay,k),a[k+'_final']) for k in STATE_KEYS)
        grid=FluidSampler(solver.mesh,record['spec'].get('degree',5),18)
        energy_errors=[];history_errors=[];norms=[]
        for i,t in enumerate(a['times']):
            s=state('_snapshots',i);metrics,energy=grid.measure(s,ref,float(t));norms.append(metrics)
            energy_errors.append(float(np.linalg.norm(energy-a['physical_energy_snapshots'][i])/np.linalg.norm(energy)))
            actual=np.array(totals(solver,s,float(t),a['spatial_ledger_snapshots'][i]))
            index=int(round(t/record['dt']));target=a['history'][index]
            # Last two fields are transient solver residuals, not state totals.
            history_errors.append(float(np.max(abs(actual[:-2]-target[:-2])/np.maximum(abs(target[:-2]),1e-30))))
        rate_norms,energy_rate=grid.rates(initial,first_rate,ref,0.)
        error=float(np.linalg.norm(energy_rate-a['initial_physical_energy_rate'])/max(np.linalg.norm(energy_rate),1e-300))
        assert max(energy_errors)<2e-11 and max(history_errors)<1e-13
        # Absolute conditioning is retained for an equilibrium's tiny residual.
        rate_rms=float(np.sqrt(np.sum(solver.mesh.volumes*(energy_rate-a['initial_physical_energy_rate'])**2)/solver.mesh.volumes.sum()))
        item=dict(file=path.name,status=record['status'],initial_and_rate_bitwise=True,final_step_bitwise=record['status']=='complete',
                  max_independent_local_energy_relative=max(energy_errors),max_history_relative=max(history_errors),
                  initial_energy_rate_relative=error,initial_energy_rate_rms_difference=rate_rms,
                  independent_order18_norms=norms,independent_order18_rate_norms=rate_norms)
        records.append(item);print(path.name,'audited',flush=True)
    result=dict(source_sha256=identity['source_sha256'],records=records)
    (directory/'audit.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);audit(parser.parse_args().directory)
