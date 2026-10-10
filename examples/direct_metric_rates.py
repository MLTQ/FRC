"""Complete signed workflows with direct construction of the used metric."""
from . import streamed_radial_rates as driver
from .cached_streamed_rates import interfaces
from .implicit_plasma_records import provenance
from plasma_bench.direct_metric_wall_plasma import direct_metric


def run(destination, specs, **options):
    original = driver.split_inverse, driver.interfaces, driver.provenance
    try:
        driver.split_inverse = direct_metric
        driver.interfaces = interfaces
        driver.provenance = lambda seeds: provenance([*seeds, 'examples.direct_metric_rates',
            'examples.streamed_interface_audit', 'examples.streamed_radial_rates'])
        return driver.run(destination, specs, **options)
    finally:
        driver.split_inverse, driver.interfaces, driver.provenance = original
