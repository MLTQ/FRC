"""Original signed physical audit with the actual split-inverse full RHS."""
from . import matched_mode_audit as audit
from . import matched_mode_comparison as comparison
from .matched_mode_fixture import fixture as original_fixture
from .implicit_plasma_records import provenance
from plasma_bench.split_inverse_wall_plasma import split_inverse


def fixture(spec):
    values = list(original_fixture(spec))
    values[3] = split_inverse(values[3])
    return tuple(values)


def run(negative, positive, destination):
    original = audit.fixture, audit.provenance
    try:
        audit.fixture = fixture
        audit.provenance = lambda seeds: provenance([*seeds, 'examples.split_mode_audit'])
        return audit.run(negative, positive, destination)
    finally:
        audit.fixture, audit.provenance = original


def compare(first, second, destination):
    original = comparison.fixture, comparison.provenance
    try:
        comparison.fixture = fixture
        comparison.provenance = lambda seeds: provenance([*seeds, 'examples.split_mode_audit'])
        return comparison.run(first, second, destination)
    finally:
        comparison.fixture, comparison.provenance = original
