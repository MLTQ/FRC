"""Supplemental stronger-gradient and common-time radial controls."""
from pathlib import Path
import hashlib
from time import perf_counter
from .warm_cases import protocol as base_protocol
from .warm_measurements import measure
from .warm_records import provenance, save
from .verification_records import write_record


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def protocol():
    base = {c['name']: c for c in base_protocol()}
    cases = []
    for kind in ('contact', 'simple_wave'):
        for nz in (128, 256):
            cases.append(base[f'{kind}_strong_64'] | dict(name=f'{kind}_strong_{nz}', nz=nz))
    for nr in (4, 8):
        cases.append(base['contact_64'] | dict(name=f'contact_common_time_nr{nr}', nr=nr, periods=.2, fixed_steps=512))
    return cases


def run(directory='outputs/warm_adiabatic_controls'):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source, control_hash, records = provenance(), source_hash(), []
    for case in protocol():
        start = perf_counter()
        record, arrays = measure(case)
        if source_hash() != control_hash:
            raise RuntimeError('Supplemental protocol source changed during run')
        records.append(save(directory, case, source, record | dict(control_source_sha256=control_hash, elapsed_seconds=perf_counter()-start), arrays))
        print(case['name'], record['metrics'], flush=True)
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, control_source_sha256=control_hash, cases=records))


if __name__ == '__main__':
    run()
