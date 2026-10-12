"""One complete physical field-refinement workflow with own-process resources."""
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import threading
import time
from .direct_metric_rates import run as calculate
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze


def run(destination, control, nr, nz, original=None):
    destination = Path(destination); control = Path(control)
    if destination.exists() or control.exists(): raise ValueError('Fresh calculation and control paths required')
    identity = provenance(['examples.guarded_field_refinement', 'examples.streamed_radial_rates',
        'examples.streamed_interface_audit', 'examples.cached_streamed_rates'])
    control.mkdir(parents=True); freeze(control, identity)
    spec = dict(elements=64, points=256, field_nr=nr, field_nz=nz, epsilon=1e-4, tolerance=1e-10)
    state = dict(terminal=False, source_sha256=identity['source_sha256'], destination=str(destination),
        original=str(original) if original is not None else None, spec=spec, pid=os.getpid(),
        disk_guard_triggered=False, minimum_free_bytes=shutil.disk_usage('.').free,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        threshold_free_bytes=256*1024**2, physical_workflow_terminal=False)
    def save(): (control/'summary.json').write_text(json.dumps(state, indent=2)+'\n')
    def observe():
        free = shutil.disk_usage('.').free
        state['minimum_free_bytes'] = min(state['minimum_free_bytes'], free)
        state['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return free
    done = threading.Event()
    def guard():
        while not done.wait(1):
            free = observe()
            if free < state['threshold_free_bytes']:
                state.update(disk_guard_triggered=True, free_bytes_at_interruption=free); save()
                os.kill(os.getpid(), signal.SIGINT); return
    save(); start = time.monotonic(); worker = threading.Thread(target=guard, daemon=True); worker.start()
    try:
        result = calculate(destination, [spec])
        state.update(physical_workflow_terminal=result['terminal'], all_cases_complete=result['all_cases_complete'],
            physical_source_sha256=result['source_sha256'])
        if original is not None:
            verify(original, destination, str(destination)+'_verification')
            state['complete_replay_verified'] = True
    finally:
        done.set(); worker.join(); observe()
        state.update(terminal=True, elapsed_seconds=time.monotonic()-start,
            rss_units='macOS own-process cumulative maximum bytes',
            scope='Resource guard changes no physical calculation or numerical acceptance; interrupted data remains incomplete')
        if provenance(['examples.guarded_field_refinement', 'examples.streamed_radial_rates',
            'examples.streamed_interface_audit', 'examples.cached_streamed_rates']) != identity:
            state['source_changed'] = True; save(); raise ValueError('Source changed')
        save()
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination'); parser.add_argument('control')
    parser.add_argument('--nr', type=int, required=True); parser.add_argument('--nz', type=int, required=True)
    parser.add_argument('--original')
    args = parser.parse_args()
    run(args.destination, args.control, args.nr, args.nz, args.original)
