"""Complete tensor workflow with identical current-work expressions/lower lifetimes."""
from . import tensor_metric_rates as tensor
from plasma_bench import regular_plasma, direct_metric_wall_plasma
from plasma_bench.released_current_transport import ReleasedCurrentTransport
from plasma_bench.released_regular_current import ReleasedRegularCurrent
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    saved = (regular_plasma.RegularCurrent, direct_metric_wall_plasma.CurrentTransport, tensor.provenance)
    try:
        regular_plasma.RegularCurrent = ReleasedRegularCurrent
        direct_metric_wall_plasma.CurrentTransport = ReleasedCurrentTransport
        tensor.provenance = lambda seeds: provenance([*seeds, 'examples.released_metric_rates',
            'examples.streamed_radial_rates', 'examples.cached_streamed_rates',
            'examples.streamed_interface_audit', 'plasma_bench.regular_plasma',
            'plasma_bench.direct_metric_wall_plasma'])
        return tensor.run(destination, specs, **options)
    finally:
        regular_plasma.RegularCurrent, direct_metric_wall_plasma.CurrentTransport, tensor.provenance = saved
