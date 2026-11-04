"""Recover and continue the original HOPS campaign with an exact replay gate."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from examples.hybrid3d_campaign import freeze_source, write_json
from examples.iec_3d import run as simulate
from examples.iec_recover import recover, exact
from plasma_bench.iec_checkpoint import load_checkpoint


def worker(args):
    out,experiment,campaign=map(Path,(args.out,args.run,args.campaign))
    lock=(campaign/'.continuation.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    progress=dict(status='recovering',pid=os.getpid(),run=str(experiment))
    write_json(out/'progress.json',progress)
    config=json.loads((experiment/'run.json').read_text())['config']
    checkpoint=experiment/'checkpoint.npz'
    if not checkpoint.exists():
        config,previous,last,paths=recover(experiment,out)
        progress['status']='replaying_saved_interval';write_json(out/'progress.json',progress)
        proof=out/'replay'
        simulate(proof,config,experiment/'mesh.npz',checkpoint=paths[previous],stop_step=last)
        old=json.loads((experiment/'history.json').read_text())[-1]
        new=json.loads((proof/'history.json').read_text())[-1]
        def physical(row): return {k:v for k,v in row.items() if k!='elapsed_s'}
        original=sorted(experiment.glob('state_*.npz'))[-1]
        replay=sorted(proof.glob('state_*.npz'))[-1]
        with np.load(original) as a,np.load(replay) as b:
            checks={k:exact(a[k],b[k]) for k in a.files}
        recovered_meta,_,_=load_checkpoint(paths[last],config,experiment/'mesh.npz')
        replay_meta,_,_=load_checkpoint(proof/'checkpoint.npz',config,experiment/'mesh.npz')
        report=dict(original=str(original),replay=str(replay),all_arrays_bitwise_identical=checks,
                    physical_diagnostics_identical=physical(old)==physical(new),
                    rng_state_identical=recovered_meta['rng']==replay_meta['rng'],
                    controller_state_identical=recovered_meta['drive']==replay_meta['drive'])
        write_json(out/'replay_verification.json',report)
        if not (all(checks.values()) and report['physical_diagnostics_identical'] and
                report['rng_state_identical'] and report['controller_state_identical']):
            raise RuntimeError('Saved-interval replay differs; original remains interrupted')
        shutil.copy2(paths[last],checkpoint)
    c=json.loads((campaign/'campaign.json').read_text())
    c.update(status='running',active_coordinator_pid=os.getpid(),continuation=str(out))
    c.pop('error',None)
    root=experiment.parents[2]
    for mode in ('hops','dc','fixed'):
        if (campaign/'STOP').exists():
            c['status']='stopped';break
        target=experiment if mode=='hops' else root/'outputs/lab_runs'/f'{campaign.name}_{mode}'
        cfg=config if mode=='hops' else config|dict(drive_mode=mode,run_label=f'HOPS · {mode} · +3 kV · 0–5 µs')
        existing=json.loads((target/'run.json').read_text()) if (target/'run.json').exists() else None
        if existing and existing['status']=='completed':
            c['cases'][mode]=dict(status='completed',run=str(target),reached_us=existing['time_ns']/1000)
            continue
        c['cases'][mode]=dict(status='running',run=str(target),target_us=cfg['duration_ns']/1000)
        write_json(campaign/'campaign.json',c)
        progress.update(status='running',mode=mode,run=str(target));write_json(out/'progress.json',progress)
        if existing:
            cp=target/'checkpoint.npz'
            if not cp.exists(): raise RuntimeError('Existing case has no complete checkpoint')
            result=simulate(target,cfg,checkpoint=cp,resume=True)
        else:
            result=simulate(target,cfg,campaign/'mesh.npz')
        c['cases'][mode].update(status=result['status'],reached_us=result['time_ns']/1000)
        write_json(campaign/'campaign.json',c)
        if result['status']!='completed':
            c['status']='stopped';break
    else:
        c['status']='completed'
    write_json(campaign/'campaign.json',c)
    progress['status']=c['status'];write_json(out/'progress.json',progress)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out',type=Path)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:
        try: return worker(args)
        except BaseException as error:
            write_json(args.out/'progress.json',dict(status='failed',pid=os.getpid(),error=str(error)))
            p=args.campaign/'campaign.json';c=json.loads(p.read_text())
            if c.get('active_coordinator_pid')==os.getpid():
                c.update(status='failed',error=str(error));write_json(p,c)
            raise
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=freeze_source(Path(__file__).resolve().parents[1],out)
    command=[sys.executable,'-m','examples.hops_continue',str(out),'--worker',
             '--run',str(args.run.resolve()),'--campaign',str(args.campaign.resolve())]
    env=os.environ|dict(PYTHONPATH=str(source),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
    with (out/'worker.log').open('w') as log:
        child=subprocess.Popen(['nohup',*command],cwd=source,env=env,stdin=subprocess.DEVNULL,
                               stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    write_json(out/'launch.json',dict(pid=child.pid,command=command,source=str(source),detached=True))
    print(json.dumps(dict(pid=child.pid,progress=str(out/'progress.json'))),flush=True)


if __name__=='__main__': main()
