"""Full released-workspace calculation with detached geometry-cache templates."""
from . import direct_metric_rates as direct
from . import released_metric_rates as released
from plasma_bench.detached_metric_wall_plasma import detached_metric
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    saved = direct.direct_metric, released.provenance
    try:
        direct.direct_metric = detached_metric
        released.provenance = lambda seeds: provenance([*seeds, 'examples.detached_metric_rates',
            'examples.released_metric_rates', 'examples.direct_metric_rates',
            'examples.tensor_metric_rates', 'plasma_bench.detached_metric_wall_plasma'])
        return released.run(destination, specs, **options)
    finally:
        direct.direct_metric, released.provenance = saved
