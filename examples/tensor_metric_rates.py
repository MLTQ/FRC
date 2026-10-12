"""Opt-in full signed workflows with separable physical electric-mass solves."""
from . import direct_metric_rates as direct
from . import regular_plasma_reference as reference
from plasma_bench.tensor_mass_geometry import TensorMassGeometry
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    original = reference.CompatibleGeometry, direct.provenance
    try:
        reference.CompatibleGeometry = TensorMassGeometry
        direct.provenance = lambda seeds: provenance([*seeds, 'examples.tensor_metric_rates',
            'plasma_bench.tensor_mass_geometry', 'examples.regular_plasma_reference'])
        return direct.run(destination, specs, **options)
    finally:
        reference.CompatibleGeometry, direct.provenance = original
