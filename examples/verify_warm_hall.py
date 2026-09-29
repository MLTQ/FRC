"""Audit all warm Hall cases, compare trajectories and render scientific evidence."""
from pathlib import Path
import hashlib
from .warm_hall_audit import audit
from .warm_hall_comparisons import comparisons
from .warm_hall_figures import convergence, controls
from .verification_records import write_record


def run(directory='outputs/warm_hall_verification'):
    directory = Path(directory)
    data, records, arrays, meshes, audits = audit(directory)
    differences = comparisons(records, arrays, meshes)
    figures = [convergence(directory, records), controls(directory, records, arrays, differences)]
    sources = ['examples/'+name+'.py' for name in ('verify_warm_hall', 'warm_hall_audit', 'warm_hall_field_audit',
               'warm_hall_comparisons', 'warm_hall_figures', 'warm_audit', 'compatible_hall_audit')]
    report = dict(completed=True, provenance=data['provenance'], audits=audits, comparisons=differences,
                  report_sources={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},
                  figures={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in figures},
                  scope='Equal-gamma massless-electron warm Hall verification; no resistive heating closure or vacuum formation')
    write_record(directory/'report.json', report)
    print(dict(cases=len(records), source=data['provenance']['source_sha256'],
               time_ratios={k: v['pressure_difference_ratios'] for k, v in differences['time'].items()}), flush=True)


if __name__ == '__main__':
    run()
