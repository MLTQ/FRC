"""Complete detached-cache signed workflow with independent-work scheduling."""
from . import detached_metric_rates as detached
from . import released_metric_rates as released
from plasma_bench.scheduled_current_transport import ScheduledCurrentTransport
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    saved = released.ReleasedCurrentTransport, detached.provenance
    try:
        released.ReleasedCurrentTransport = ScheduledCurrentTransport
        detached.provenance = lambda seeds: provenance([*seeds, 'examples.scheduled_metric_rates',
            'examples.detached_metric_rates', 'examples.released_metric_rates',
            'plasma_bench.scheduled_current_transport'])
        return detached.run(destination, specs, **options)
    finally:
        released.ReleasedCurrentTransport, detached.provenance = saved
