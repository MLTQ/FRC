"""Guard the complete three-grid coupled response comparison and replay."""
import argparse
import importlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import threading
import time
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze


def run(kind, destination, control, inputs=(), original=None):
    names={'comparison':'examples.coupled_refinement_comparison'}
    if kind not in names: raise ValueError('Known spatial contribution workflow required')
    if len(inputs)!=0: raise ValueError('Exact declared inputs required')
    destination=Path(destination);control=Path(control)
    if destination.exists() or control.exists(): raise ValueError('Fresh output/control paths required')
    seeds=['examples.guarded_coupled_comparison',names[kind]]
    identity=provenance(seeds);control.mkdir(parents=True);freeze(control,identity)
    state=dict(terminal=False,pid=os.getpid(),kind=kind,destination=str(destination),inputs=list(inputs),
               original=original,source_sha256=identity['source_sha256'],workflow_terminal=False,
               resource_guard_triggered=False,minimum_free_bytes=shutil.disk_usage('.').free,peak_rss_bytes=0,
               threshold_free_bytes=256*1024**2,maximum_rss_bytes=8*1024**3)
    def save():
        p=control/'summary.pending.json';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(control/'summary.json')
    def observe():
        state['minimum_free_bytes']=min(state['minimum_free_bytes'],shutil.disk_usage('.').free)
        state['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return state['minimum_free_bytes']>=state['threshold_free_bytes'] and state['peak_rss_bytes']<=state['maximum_rss_bytes']
    done=threading.Event()
    def guard():
        count=0
        while not done.wait(1):
            count+=1
            if not observe():
                state['resource_guard_triggered']=True;save();os.kill(os.getpid(),signal.SIGINT)
                if not done.wait(5):
                    state.update(terminal=True,forced_own_process_exit=True,interruption='Own worker exits130 after native work did not yield within5s')
                    save();os._exit(130)
                return
            if count%30==0:
                state.update(monitor_observations=count,observed_elapsed_seconds=time.monotonic()-started);save()
    if not observe():
        state.update(terminal=True,resource_guard_triggered=True,interruption='Resource preflight rejected');save()
        raise RuntimeError('Insufficient resources before coupled comparison')
    save();started=time.monotonic();worker=threading.Thread(target=guard,daemon=True);worker.start()
    try:
        result=importlib.import_module(names[kind]).run(destination,*inputs)
        state.update(workflow_terminal=result['terminal'],all_cases_complete=result.get('all_cases_complete',True),
                     scientific_source_sha256=result['source_sha256'])
        if original is not None:
            target=str(destination)+'_verification';verify(original,destination,target)
            replay=json.loads((Path(target)/'summary.json').read_text())
            if any(v['original_sha256']!=v['repeated_sha256'] for v in replay['arrays'].values()):
                raise ValueError('Every raw file must replay byte identically')
            state['complete_protocol_replay_verified']=True
    except BaseException as error:
        state.update(exception_type=type(error).__name__,exception_message=str(error));raise
    finally:
        done.set();worker.join();within=observe()
        state.update(terminal=True,resources_within_bounds=within,elapsed_seconds=time.monotonic()-started,
                     rss_units='macOS own-process cumulative maximum bytes')
        if provenance(seeds)!=identity:
            state['source_changed']=True;save();raise ValueError('Coupled comparison source changed')
        save()
    if not within or state['resource_guard_triggered']:raise RuntimeError('Recorded resource limits exceeded')
    return state


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('kind',choices=('comparison',));p.add_argument('destination')
    p.add_argument('control');p.add_argument('inputs',nargs='*');p.add_argument('--original')
    a=p.parse_args();run(a.kind,a.destination,a.control,a.inputs,a.original)
