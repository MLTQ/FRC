"""Serial isolated full trajectories with resource-budgeted fresh replays."""
import argparse
import json
from pathlib import Path
import os
import shutil
import signal
import subprocess
import sys
from .refined_control_campaign import share_verified_replay
from .implicit_plasma_records import provenance, freeze, sha


def storage_reserve(spec, steps):
    if (set(spec)!={'elements','points','field_nr','field_nz'}
            or any(type(v) is not int or v<1 for v in spec.values())
            or type(steps) is not int or steps<1):
        raise ValueError('Four positive integer mesh counts and positive step count required')
    e,p,nr,nz=(spec[k] for k in ('elements','points','field_nr','field_nz'))
    if p!=nz:raise ValueError('Storage derivation requires the declared matching axial counts')
    state=36*e*p;magnetic=(2*nr+1)*nz
    coordinates=(36*e-3)*p+2*(nr-1)*nz
    initial_bytes=8*(state+magnetic)
    step_bytes=initial_bytes+8*2*coordinates
    payload=3*initial_bytes+steps*step_bytes
    return dict(uncompressed_payload_bytes=payload,raw_step_bytes=step_bytes,
        reserved_bytes=(11*payload+9)//10+8*1024**2,
        scope='1.1x uncompressed payload plus8MiB metadata/archive planning reserve; not a guaranteed bound or a Git estimate')


def run(destination, plan_path):
    destination=Path(destination);plan_path=Path(plan_path)
    if destination.exists():raise ValueError('Fresh campaign path required')
    jobs=json.loads(plan_path.read_text());plan_hash=sha(plan_path)
    seeds=['examples.trajectory_control_campaign','examples.guarded_trajectory_controls',
        'examples.split_mode_evolution','examples.matched_mode_evolution','examples.matched_mode_fixture']
    identity=provenance(seeds);destination.mkdir(parents=True);freeze(destination,identity)
    state=dict(terminal=False,pid=os.getpid(),source_sha256=identity['source_sha256'],
        plan=str(plan_path),plan_sha256=plan_hash,jobs=[],shared_verified_files=[])
    def save():
        temporary=destination/'summary.pending.json'
        temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(destination/'summary.json')
    save()
    try:
        for job in jobs:
            if not job['name'].isidentifier():raise ValueError('Simple job name required')
            budget=storage_reserve(job['spec'],job['steps']);free=shutil.disk_usage('.').free
            row=dict(job=job,status='preflight',storage_budget=budget,free_bytes_before=free)
            state['jobs'].append(row);save()
            if free<256*1024**2+budget['reserved_bytes']:
                row['status']='not_started_resource_preflight';save()
                raise RuntimeError('Insufficient estimated disk reserve before '+job['name'])
            if job.get('original'):
                before=json.loads((Path(job['original'])/'summary.json').read_text())
                if not before['terminal']:raise ValueError('Terminal original protocol required before replay')
            args=[sys.executable,'-m','examples.guarded_trajectory_controls',job['destination'],job['control'],
                '--spec',json.dumps(job['spec'],sort_keys=True),'--epsilon='+str(job['epsilon']),
                '--steps',str(job['steps'])]
            if job.get('original'):args+=['--original',job['original']]
            log=destination/(job['name']+'.log');row.update(log=str(log),status='starting');save()
            print('Starting',job['name'],flush=True)
            with log.open('xb') as output:
                child=subprocess.Popen(args,stdout=output,stderr=subprocess.STDOUT)
                row.update(pid=child.pid,status='running');save()
                try:code=child.wait()
                except BaseException:
                    if child.poll() is None:
                        child.send_signal(signal.SIGINT)
                        try:child.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            child.terminate();child.wait(timeout=10)
                    row.update(returncode=child.returncode,status='interrupted',log_sha256=sha(log));save()
                    raise
            row.update(returncode=code,status='returned',log_sha256=sha(log));save()
            control=Path(job['control'])/'summary.json'
            if control.exists():row.update(control_sha256=sha(control),control_record=json.loads(control.read_text()))
            if code!=0:save();raise RuntimeError('Child interrupted or failed: '+job['name'])
            c=row['control_record']
            if not c['terminal'] or not c['physical_workflow_terminal'] or not c['resources_within_bounds'] or c['resource_guard_triggered']:
                raise ValueError('Terminal resource-compliant trajectory protocol required')
            row['status']='full_duration_complete' if c['full_duration_complete'] else 'trajectory_failure'
            if job.get('original'):
                if not c.get('complete_protocol_replay_verified'):raise ValueError('Complete protocol replay required')
                state['shared_verified_files'].extend(share_verified_replay(job['destination']))
            save();print('Finished',job['name'],row['status'],flush=True)
        state.update(all_jobs_terminal=True,
            all_full_duration_complete=all(r['control_record']['full_duration_complete'] for r in state['jobs']))
    finally:
        state['terminal']=True
        if sha(plan_path)!=plan_hash or provenance(seeds)!=identity:
            state['source_or_plan_changed']=True;save();raise ValueError('Source or plan changed')
        save()
    return state


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('plan')
    args=parser.parse_args();run(args.destination,args.plan)
