"""Resource-bounded tensor candidate with escalation inside its own worker."""
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import threading
import time
from .tensor_metric_rates import run as calculate
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze


def run(destination, control, nr, nz, original=None):
    destination = Path(destination); control = Path(control)
    if destination.exists() or control.exists(): raise ValueError('Fresh output/control paths required')
    seeds = ['examples.guarded_tensor_refinement', 'examples.streamed_radial_rates',
        'examples.cached_streamed_rates', 'examples.streamed_interface_audit']
    identity = provenance(seeds); control.mkdir(parents=True); freeze(control, identity)
    spec = dict(elements=64, points=256, field_nr=nr, field_nz=nz, epsilon=1e-4, tolerance=1e-10)
    state = dict(terminal=False, source_sha256=identity['source_sha256'], destination=str(destination),
        original=original, spec=spec, pid=os.getpid(), resource_guard_triggered=False,
        minimum_free_bytes=shutil.disk_usage('.').free, peak_rss_bytes=0,
        threshold_free_bytes=256*1024**2, maximum_rss_bytes=8*1024**3, physical_workflow_terminal=False)
    def save(): (control/'summary.json').write_text(json.dumps(state, indent=2)+'\n')
    def observe():
        free = shutil.disk_usage('.').free; rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        state['minimum_free_bytes'] = min(state['minimum_free_bytes'], free); state['peak_rss_bytes'] = rss
        return free, rss
    def violation():
        free, rss = observe()
        return free < state['threshold_free_bytes'] or rss > state['maximum_rss_bytes']
    done = threading.Event()
    def guard():
        while not done.wait(1):
            if violation():
                state['resource_guard_triggered'] = True; save()
                os.kill(os.getpid(), signal.SIGINT)
                if not done.wait(5):
                    state.update(terminal=True, forced_own_process_exit=True,
                        interruption='Native work did not yield to SIGINT within5s; own worker exits130')
                    save(); os._exit(130)
                return
    observe(); save()
    if violation():
        state.update(terminal=True, resource_guard_triggered=True, interruption='Insufficient headroom before initialization')
        save(); raise RuntimeError('Resource preflight rejected before physical work')
    start = time.monotonic(); worker = threading.Thread(target=guard, daemon=True); worker.start()
    try:
        result = calculate(destination, [spec])
        state.update(physical_workflow_terminal=result['terminal'], all_cases_complete=result['all_cases_complete'],
            physical_source_sha256=result['source_sha256'])
        if original is not None:
            verify(original, destination, str(destination)+'_verification'); state['complete_replay_verified'] = True
    finally:
        done.set(); worker.join(); observe()
        state.update(terminal=True, elapsed_seconds=time.monotonic()-start,
            rss_units='macOS own-process cumulative maximum bytes',
            scope='Resource interruption is not numerical rejection; no altered physical acceptance')
        if provenance(seeds) != identity:
            state['source_changed'] = True; save(); raise ValueError('Source changed')
        save()
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination'); parser.add_argument('control')
    parser.add_argument('--nr', type=int, required=True); parser.add_argument('--nz', type=int, required=True)
    parser.add_argument('--original'); args = parser.parse_args()
    run(args.destination, args.control, args.nr, args.nz, args.original)
