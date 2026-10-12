"""Complete scheduled signed workflow with componentwise drift stress."""
from . import scheduled_metric_rates as scheduled
from plasma_bench.component_stress_transport import ComponentStressTransport
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    saved = scheduled.ScheduledCurrentTransport, scheduled.provenance
    try:
        scheduled.ScheduledCurrentTransport = ComponentStressTransport
        scheduled.provenance = lambda seeds: provenance([*seeds, 'examples.component_stress_rates',
            'examples.scheduled_metric_rates', 'plasma_bench.component_stress_transport'])
        return scheduled.run(destination, specs, **options)
    finally:
        scheduled.ScheduledCurrentTransport, scheduled.provenance = saved
