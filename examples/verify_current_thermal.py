"""Read-only numerical audit and figures for the complete current/thermal study."""
import argparse
from pathlib import Path
import hashlib
import json
from .current_thermal_audit import collection
from .current_thermal_comparisons import comparisons
from .current_thermal_figures import figures
from .current_symbol_measurement import measure as symbol_measure
from .verification_records import write_record


def verify(directory):
    records,archives,meshes,audits=collection(directory)
    print('Audited all',len(records),'archives and accepted histories',flush=True)
    compared=comparisons(records,archives,meshes)
    symbol=json.loads((directory/'principal_symbol.json').read_text())
    assert {k:v for k,v in symbol.items() if k!='provenance'}==symbol_measure()
    paths=figures(directory,records,archives,meshes,compared,symbol)
    source_names=('examples/verify_current_thermal.py','examples/current_thermal_audit.py','examples/current_thermal_comparisons.py',
        'examples/current_thermal_figures.py','tests/test_current_thermal.py','tests/test_toroidal_thermal_symbol.py')
    summary=dict(total_records=len(records),succeeded=sum(r['completed'] for r in records.values()),
        rejected=[k for k,r in records.items() if not r['completed']],
        max_component_energy_relative_drift=max(a['max_total_energy_j_relative_drift'] for a in audits.values()),
        max_conservative_energy_relative_drift=max(a['max_conservative_energy_j_relative_drift'] for a in audits.values()),
        max_electron_log_entropy_change_per_particle=max(a['electron_max_log_entropy_change_per_particle'] for a in audits.values()),
        max_ion_log_entropy_change_per_particle=max(a['ion_max_log_entropy_change_per_particle'] for a in audits.values()))
    report=dict(audit_passed=True,engine_provenance=next(iter(records.values()))['provenance'],summary=summary,
        report_sources={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in source_names},
        figures={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},audits=audits,comparisons=compared,
        principal_symbol=symbol,scope='Physical finite-mass, fixed n and constrained zero barycentric flow; no moving formation validation.')
    write_record(directory/'report.json',report);print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/current_thermal_verification')
    verify(Path(parser.parse_args().directory))
