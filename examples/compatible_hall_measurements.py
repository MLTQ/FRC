"""Run the declared compatible Hall spatial and actual fluid-evolution protocol."""
from pathlib import Path
from .compatible_hall_cases import protocol
from .compatible_hall_records import provenance
from .compatible_hall_spatial import measure as spatial
from .compatible_hall_evolution import measure as evolve
from .verification_records import write_record


def run(directory='outputs/compatible_hall_verification'):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    records = [(spatial if case['kind'] == 'spatial' else evolve)(case, directory) for case in protocol()]
    if provenance() != source:
        raise RuntimeError('Compatible Hall source changed across protocol')
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, cases=records))


if __name__ == '__main__':
    run()
