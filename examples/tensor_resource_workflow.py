"""Read-only stage resource observations during a complete fine baseline replay."""
import json
from pathlib import Path
import resource
import time
from plasma_bench.current_transport import CurrentTransport
from plasma_bench.physical_derivatives import CurrentDifferential, RegularDifferential
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.direct_metric_wall_plasma import DirectMetricWallInduction
from . import guarded_tensor_refinement as guarded
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, observed):
    destination = Path(destination); observed = Path(observed)
    if destination.exists() or observed.exists(): raise ValueError('Fresh paths required')
    seeds = ['examples.tensor_resource_workflow', 'examples.guarded_tensor_refinement',
             'examples.streamed_radial_rates', 'examples.cached_streamed_rates',
             'examples.streamed_interface_audit']
    identity = provenance(seeds); observed.mkdir(parents=True); freeze(observed, identity)
    calls = []; patches = []
    targets = ((CurrentTransport, 'pair'), (CurrentTransport, 'advection'),
               (CurrentDifferential, 'gradient'), (RegularDifferential, 'velocity_gradient'),
               (RegularDifferential, 'stress_rate'), (RegularCurrent, 'electromagnetic_pair'),
               (RegularSpecies, 'pair'), (DirectMetricWallInduction, '__init__'))
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
    for cls, name in targets:
        old = getattr(cls, name); patches.append((cls, name, old))
        setattr(cls, name, wrap(old, cls.__name__+'.'+name))
    try:
        control = guarded.run(destination, observed/'guard', 256, 256,
                              original='outputs/tensor_mass_fine_baseline')
    finally:
        for cls, name, method in reversed(patches): setattr(cls, name, method)
    if not control.get('complete_replay_verified'): raise ValueError('Complete scientific replay required')
    if provenance(seeds) != identity: raise ValueError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], destination=str(destination),
                  calls_file='calls.json', calls_sha256=sha(observed/'calls.json'),
                  full_replay_verified=True, resource_guard_triggered=control['resource_guard_triggered'],
                  peak_rss_bytes=control['peak_rss_bytes'],
                  scope='Read-only observations. macOS ru_maxrss is cumulative process high-water memory, not live allocation or per-function exclusive storage; nested intervals overlap. Runtime/resource values need not replay. The entire signed physical calculation and every coefficient/audit must reproduce the baseline.')
    (observed/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
