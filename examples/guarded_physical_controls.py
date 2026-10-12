"""Explicit amplitude/residual controls with the verified physical workflow."""
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
from .component_stress_rates import run as calculate
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze


def run(destination, control, nr, nz, epsilon, tolerance, original=None):
    destination = Path(destination); control = Path(control)
    if destination.exists() or control.exists(): raise ValueError('Fresh physical/control paths required')
    if (type(nr) is not int or nr < 64 or nr % 64 or type(nz) is not int or nz < 6
            or not math.isfinite(epsilon) or epsilon <= 0
            or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-10):
        raise ValueError('Resolved F64 field subdivision, positive amplitude and unchanged/tighter bound required')
    seeds = ['examples.guarded_physical_controls', 'examples.component_stress_rates',
        'examples.scheduled_metric_rates', 'examples.detached_metric_rates',
        'examples.released_metric_rates', 'examples.tensor_metric_rates',
        'examples.direct_metric_rates', 'examples.streamed_radial_rates',
        'examples.cached_streamed_rates', 'examples.streamed_interface_audit',
        'examples.regular_plasma_reference']
    identity = provenance(seeds); control.mkdir(parents=True); freeze(control, identity)
    spec = dict(elements=64, points=256, field_nr=nr, field_nz=nz, epsilon=epsilon, tolerance=tolerance)
    state = dict(terminal=False, source_sha256=identity['source_sha256'], destination=str(destination),
        original=original, spec=spec, pid=os.getpid(), resource_guard_triggered=False,
        minimum_free_bytes=shutil.disk_usage('.').free, peak_rss_bytes=0,
        threshold_free_bytes=256*1024**2, maximum_rss_bytes=8*1024**3, physical_workflow_terminal=False)
    def save(): (control/'summary.json').write_text(json.dumps(state, indent=2)+'\n')
    def observe():
        state['minimum_free_bytes'] = min(state['minimum_free_bytes'], shutil.disk_usage('.').free)
        state['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return (state['minimum_free_bytes'] >= state['threshold_free_bytes']
                and state['peak_rss_bytes'] <= state['maximum_rss_bytes'])
    done = threading.Event()
    def guard():
        while not done.wait(1):
            if not observe():
                state['resource_guard_triggered'] = True; save(); os.kill(os.getpid(), signal.SIGINT)
                if not done.wait(5):
                    state.update(terminal=True, forced_own_process_exit=True,
                        interruption='Native work did not yield to SIGINT within5s; own worker exits130')
                    save(); os._exit(130)
                return
    if not observe():
        state.update(terminal=True, resource_guard_triggered=True, interruption='Resource preflight rejected')
        save(); raise RuntimeError('Insufficient resources before physical initialization')
    save(); started = time.monotonic(); worker = threading.Thread(target=guard, daemon=True); worker.start()
    try:
        result = calculate(destination, [spec])
        state.update(physical_workflow_terminal=result['terminal'], all_cases_complete=result['all_cases_complete'],
            physical_source_sha256=result['source_sha256'],
            numerical_rejection=any(case['status'] == 'failed' for case in result['cases']))
        if original is not None:
            verify(original, destination, str(destination)+'_verification')
            state.update(complete_protocol_replay_verified=True,
                accepted_physical_replay=result['all_cases_complete'])
    finally:
        done.set(); worker.join(); within = observe()
        state.update(terminal=True, resources_within_bounds=within, elapsed_seconds=time.monotonic()-started,
            rss_units='macOS own-process cumulative maximum bytes',
            scope='Accepted full physical workflow, strict numerical rejection and resource interruption are distinct. A replayed rejection is not an accepted physical rate.')
        if provenance(seeds) != identity:
            state['source_changed'] = True; save(); raise ValueError('Source changed')
        save()
    if not within or state['resource_guard_triggered']: raise RuntimeError('Recorded resource limits exceeded')
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('control')
    parser.add_argument('--nr', type=int, required=True); parser.add_argument('--nz', type=int, required=True)
    parser.add_argument('--epsilon', type=float, required=True); parser.add_argument('--tolerance', type=float, required=True)
    parser.add_argument('--original'); args = parser.parse_args()
    run(args.destination, args.control, args.nr, args.nz, args.epsilon, args.tolerance, args.original)
