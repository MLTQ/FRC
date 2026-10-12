"""Observe stage high-water RSS during a complete scheduled fine replay."""
import argparse
import json
from pathlib import Path
import resource
import time
from plasma_bench.scheduled_current_transport import ScheduledCurrentTransport
from plasma_bench.released_current_transport import ReleasedCurrentTransport
from plasma_bench.released_regular_current import ReleasedRegularCurrent
from plasma_bench.physical_derivatives import CurrentDifferential, RegularDifferential
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.direct_metric_wall_plasma import DirectMetricWallInduction
from plasma_bench.split_inverse_wall_induction import SplitInverseWallInduction
from . import guarded_scheduled_refinement as guarded
from . import streamed_radial_rates as driver
from . import direct_metric_rates as direct
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, observed):
    destination = Path(destination); observed = Path(observed)
    if destination.exists() or observed.exists(): raise ValueError('Fresh paths required')
    seeds = ['examples.scheduled_resource_workflow', 'examples.guarded_scheduled_refinement',
             'examples.streamed_radial_rates', 'examples.direct_metric_rates',
             'examples.cached_streamed_rates', 'examples.streamed_interface_audit']
    identity = provenance(seeds); observed.mkdir(parents=True); freeze(observed, identity)
    calls = []; patches = []
    targets = ((ScheduledCurrentTransport, 'pair'), (ReleasedCurrentTransport, 'advection'),
               (CurrentDifferential, 'gradient'), (RegularDifferential, 'velocity_gradient'),
               (RegularDifferential, 'stress_rate'), (ReleasedRegularCurrent, 'electromagnetic_pair'),
               (RegularSpecies, 'pair'), (DirectMetricWallInduction, '__init__'),
               (SplitInverseWallInduction, 'solve'), (driver, 'measure'), (direct, 'interfaces'))
    def wrap(method, name):
        def measured(*args, **kwargs):
            record = dict(name=name, entry_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            calls.append(record); started = time.monotonic()
            try:
                result = method(*args, **kwargs); record['returned'] = True; return result
            finally:
                record.update(exit_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                              elapsed_seconds=time.monotonic()-started)
                (observed/'calls.json').write_text(json.dumps(calls, indent=2)+'\n')
        return measured
    for owner, name in targets:
        old = getattr(owner, name); owned = name in vars(owner)
        patches.append((owner, name, old, owned))
        setattr(owner, name, wrap(old, owner.__name__+'.'+name))
    try:
        control = guarded.run(destination, observed/'guard', 256, 256,
                              original='outputs/scheduled_current_fine')
    finally:
        for owner, name, method, owned in reversed(patches):
            if owned: setattr(owner, name, method)
            else: delattr(owner, name)
    if not control.get('complete_replay_verified'): raise ValueError('Complete scientific replay required')
    if (control['resource_guard_triggered'] or control['peak_rss_bytes'] > control['maximum_rss_bytes']
            or control['minimum_free_bytes'] < control['threshold_free_bytes']):
        raise ValueError('Resource bounds exceeded')
    if provenance(seeds) != identity: raise ValueError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], destination=str(destination),
        calls_file='calls.json', calls_sha256=sha(observed/'calls.json'), full_replay_verified=True,
        resource_guard_triggered=control['resource_guard_triggered'], peak_rss_bytes=control['peak_rss_bytes'],
        scope='Read-only cumulative own-process high-water RSS observations. Nested intervals overlap; these are not exclusive or live allocation measurements. All signed scientific arrays and metrics replay the unobserved baseline.')
    (observed/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('observed')
    args = parser.parse_args(); run(args.destination, args.observed)
