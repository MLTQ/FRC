"""Guard and fully replay the unchanged complete nonlinear trajectory protocol."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import threading
import time
from .split_mode_evolution import run as evolve
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze


def run(destination, control, spec, epsilon, steps, original=None):
    destination=Path(destination);control=Path(control)
    if destination.exists() or control.exists():raise ValueError('Fresh trajectory/control paths required')
    if not math.isfinite(epsilon) or epsilon==0 or type(steps) is not int or steps<1:
        raise ValueError('Finite signed amplitude and positive step count required')
    if set(spec) != {'elements','points','field_nr','field_nz'}:
        raise ValueError('Explicit four-count physical mesh required; original fixture settings remain')
    if any(type(v) is not int or v<1 for v in spec.values()):raise ValueError('Positive integer mesh counts required')
    seeds=['examples.guarded_trajectory_controls','examples.split_mode_evolution',
        'examples.matched_mode_evolution','examples.matched_mode_fixture']
    identity=provenance(seeds);control.mkdir(parents=True);freeze(control,identity)
    state=dict(terminal=False,source_sha256=identity['source_sha256'],destination=str(destination),
        original=original,spec=spec,epsilon=epsilon,steps=steps,duration=1e-4,
        absolute_tolerance=1e-12,expected_relative_tolerance=1e-10,pid=os.getpid(),
        resource_guard_triggered=False,minimum_free_bytes=shutil.disk_usage('.').free,
        peak_rss_bytes=0,threshold_free_bytes=256*1024**2,maximum_rss_bytes=8*1024**3,
        physical_workflow_terminal=False)
    def save():
        temporary=control/'summary.pending.json'
        temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(control/'summary.json')
    def observe():
        state['minimum_free_bytes']=min(state['minimum_free_bytes'],shutil.disk_usage('.').free)
        state['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return state['minimum_free_bytes']>=state['threshold_free_bytes'] and state['peak_rss_bytes']<=state['maximum_rss_bytes']
    done=threading.Event()
    def guard():
        observations=0
        while not done.wait(1):
            observations+=1
            if not observe():
                state['resource_guard_triggered']=True;save();os.kill(os.getpid(),signal.SIGINT)
                if not done.wait(5):
                    state.update(terminal=True,forced_own_process_exit=True,
                        interruption='Native work did not yield within5s; own worker exits130')
                    save();os._exit(130)
                return
            if observations%30==0:
                state.update(monitor_observations=observations,observed_elapsed_seconds=time.monotonic()-started)
                save()
    if not observe():
        state.update(terminal=True,resource_guard_triggered=True,interruption='Resource preflight rejected')
        save();raise RuntimeError('Insufficient resources before trajectory initialization')
    save();started=time.monotonic();worker=threading.Thread(target=guard,daemon=True);worker.start()
    try:
        result=evolve(destination,spec,epsilon=epsilon,steps=steps,duration=1e-4,absolute_tolerance=1e-12)
        if result['relative_tolerance'] != state['expected_relative_tolerance']:
            raise ValueError('Original relative tolerance changed')
        if (result['spec']!=spec or result['epsilon']!=epsilon or result['steps']!=steps
                or result['duration']!=state['duration'] or result['absolute_tolerance']!=state['absolute_tolerance']
                or result['final_time']!=result['completed_steps']*result['dt']):
            raise ValueError('Declared physical trajectory settings or progress changed')
        complete=result['success'] and result['completed_steps']==steps
        state.update(physical_workflow_terminal=result['terminal'],full_duration_complete=complete,
            completed_steps=result['completed_steps'],final_time=result['final_time'],
            trajectory_failure=not result['success'],failure_type=result.get('failure_type'),
            failure_message=result.get('failure_message'),physical_source_sha256=result['source_sha256'])
        if original is not None:
            path=str(destination)+'_verification';verify(original,destination,path)
            replay=json.loads((Path(path)/'summary.json').read_text())
            if any(r['original_sha256']!=r['repeated_sha256'] for r in replay['arrays'].values()):
                raise ValueError('Every raw state/stage file must be byte identical')
            state.update(complete_protocol_replay_verified=True,accepted_full_trajectory_replay=complete)
    finally:
        done.set();worker.join();within=observe()
        state.update(terminal=True,resources_within_bounds=within,elapsed_seconds=time.monotonic()-started,
            rss_units='macOS own-process cumulative maximum bytes',
            scope='Full-duration success, terminal trajectory failure and resource interruption are distinct. A replayed failure is not an accepted complete trajectory.')
        if provenance(seeds) != identity:
            state['source_changed']=True;save();raise ValueError('Source changed')
        save()
    if not within or state['resource_guard_triggered']:raise RuntimeError('Recorded resource limits exceeded')
    return state


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('control')
    parser.add_argument('--spec',required=True);parser.add_argument('--epsilon',type=float,required=True)
    parser.add_argument('--steps',type=int,required=True);parser.add_argument('--original');args=parser.parse_args()
    run(args.destination,args.control,json.loads(args.spec),args.epsilon,args.steps,args.original)
