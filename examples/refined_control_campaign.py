"""Serial isolated physical controls, with explicit child terminal records."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from .implicit_plasma_records import provenance, freeze, sha


def share_verified_replay(destination):
    """Share only complete, byte-identical, independently verified replay files."""
    verification = Path(str(destination)+'_verification')/'summary.json'
    record = json.loads(verification.read_text()); rows = []
    if not record['terminal'] or not record['bitwise_metrics']: raise ValueError('Complete replay required')
    for prefix in ('original', 'repeated'):
        if sha(Path(record[prefix+'_directory'])/'summary.json') != record[prefix+'_summary_sha256']:
            raise ValueError('Replay summary changed')
    for name, entry in record['arrays'].items():
        original = Path(record['original_directory'])/name
        repeated = Path(record['repeated_directory'])/name
        if not entry['bitwise_equal']: raise ValueError('Bitwise replay required')
        if sha(original) != entry['original_sha256'] or sha(repeated) != entry['repeated_sha256']:
            raise ValueError('Replay raw file changed')
        if entry['original_sha256'] != entry['repeated_sha256']: continue
        if original.stat().st_ino == repeated.stat().st_ino: continue
        size = repeated.stat().st_size; temporary = repeated.with_name(repeated.name+'.immutable-share')
        if temporary.exists(): raise ValueError('Unexpected sharing temporary')
        os.link(original, temporary); os.replace(temporary, repeated)
        if sha(original) != entry['original_sha256'] or sha(repeated) != entry['repeated_sha256']:
            raise ValueError('Shared content changed')
        rows.append(dict(original=str(original), repeated=str(repeated), sha256=entry['original_sha256'],
            bytes=size, verification=str(verification), verification_sha256=sha(verification),
            all_paths_and_content_preserved=True, shared_paths_immutable=True))
    return rows


def run(destination, plan_path):
    destination = Path(destination); plan_path = Path(plan_path)
    if destination.exists(): raise ValueError('Fresh campaign path required')
    jobs = json.loads(plan_path.read_text()); plan_hash = sha(plan_path)
    seeds = ['examples.refined_control_campaign', 'examples.guarded_physical_controls',
        'examples.component_stress_rates', 'examples.scheduled_metric_rates',
        'examples.detached_metric_rates', 'examples.released_metric_rates',
        'examples.tensor_metric_rates', 'examples.direct_metric_rates',
        'examples.streamed_radial_rates', 'examples.cached_streamed_rates',
        'examples.streamed_interface_audit', 'examples.regular_plasma_reference']
    identity = provenance(seeds); destination.mkdir(parents=True); freeze(destination, identity)
    state = dict(terminal=False, source_sha256=identity['source_sha256'], plan=str(plan_path),
        plan_sha256=plan_hash, pid=os.getpid(), jobs=[], shared_verified_files=[])
    def save(): (destination/'summary.json').write_text(json.dumps(state, indent=2)+'\n')
    save()
    try:
        for job in jobs:
            if not job['name'].isidentifier(): raise ValueError('Simple job name required')
            args = [sys.executable, '-m', 'examples.guarded_physical_controls', job['destination'], job['control']]
            args += [value for name in ('nr','nz','epsilon','tolerance') for value in ('--'+name,str(job[name]))]
            if job.get('original'): args += ['--original',job['original']]
            log = destination/(job['name']+'.log')
            row = dict(job=job, status='starting', log=str(log)); state['jobs'].append(row); save()
            print('Starting',job['name'],flush=True)
            with log.open('xb') as output:
                child = subprocess.Popen(args, stdout=output, stderr=subprocess.STDOUT)
                row.update(pid=child.pid, status='running'); save()
                try: code = child.wait()
                except BaseException:
                    if child.poll() is None:
                        child.send_signal(signal.SIGINT)
                        try: child.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            child.terminate(); child.wait(timeout=10)
                    row.update(returncode=child.returncode, status='interrupted', log_sha256=sha(log)); save()
                    raise
            row.update(returncode=code, status='returned', log_sha256=sha(log)); save()
            control = Path(job['control'])/'summary.json'
            if control.exists(): row.update(control_sha256=sha(control), control_record=json.loads(control.read_text()))
            if code != 0: save(); raise RuntimeError('Child interrupted or failed: '+job['name'])
            c = row['control_record']
            if not c['terminal'] or not c['physical_workflow_terminal'] or not c['resources_within_bounds'] or c['resource_guard_triggered']:
                raise ValueError('Terminal resource-compliant physical protocol required')
            row['status'] = 'accepted' if c['all_cases_complete'] else 'numerical_rejection'
            if job.get('original'):
                if not c.get('complete_protocol_replay_verified'): raise ValueError('Complete protocol replay required')
                state['shared_verified_files'].extend(share_verified_replay(job['destination']))
            save(); print('Finished',job['name'],row['status'],flush=True)
        state['all_jobs_terminal'] = True
    finally:
        state['terminal'] = True
        if sha(plan_path) != plan_hash or provenance(seeds) != identity:
            state['source_or_plan_changed'] = True; save(); raise ValueError('Source or plan changed')
        save()
    return state


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('plan')
    args=parser.parse_args();run(args.destination,args.plan)
