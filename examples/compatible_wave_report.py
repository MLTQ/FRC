"""Audit both representation studies and compare to preserved linear diagnosis."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import tarfile
from .compatible_wave_audit import source,case,digest
from .compatible_wave_verification import protocol
from .compatible_wave_plots import plot
from .verification_records import write_record


def helpers():
    seen=set()
    def visit(module):
        path=Path(module.replace('.','/')+'.py')
        if module in seen or not path.is_file():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.module:names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[a.name for a in node.names]
            else:continue
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for name in ('plasma_bench.__init__','examples.compatible_wave_report','tests.test_compatible_geometry','tests.test_harmonic_magnetic'):visit(name)
    paths=[Path(m.replace('.','/')+'.py') for m in seen]
    paths+=[p.with_suffix('.md') for p in paths if p.with_suffix('.md').is_file()]
    return {str(p):digest(p) for p in sorted(paths)}


def run(directory):
    directory=Path(directory);raw=Path('outputs/compatible_wave_rounded_attempt')
    frozen,_=source(directory,current=True);old,changed=source(raw,current=False)
    new_records=[];old_records=[];audits={};old_audits={}
    for configuration in protocol():
        name=configuration['name'];row,audit=case(directory,name,frozen)
        assert row['configuration']==configuration
        before,old_audit=case(raw,name,old)
        assert before['configuration']=={k:v for k,v in configuration.items() if k!='representation'}
        new_records.append(row);old_records.append(before);audits[name]=audit;old_audits[name]=old_audit
        print('audited affine and rounded',name,flush=True)
    previous_path=Path('outputs/smooth_finite_mass_verification/slow_mode_diagnosis.json')
    previous=json.loads(previous_path.read_text())
    assert digest(previous_path.parent/previous['fields']['file'])==previous['fields']['sha256']
    compact=lambda row:{k:v for k,v in row.items() if k not in ('provenance','elapsed_s')}
    summary=dict(completed_affine=len(new_records),completed_rounded=len(old_records),
        maximum_original_backward=max(v for a in audits.values() for v in a['backward'].values()),
        maximum_physical_energy_difference=max(a[k] for a in audits.values() for k in ('magnetic_energy_relative_difference','current_energy_relative_difference','increment_energy_relative_difference')),
        maximum_norm_quadrature_change=max(r['quadrature_absolute_change'] for r in new_records),
        maximum_full_faraday_rate_relative=max(a['full_faraday_rate_relative'] for a in audits.values()))
    report=dict(scope='Higher-order compatible field and harmonic-storage rate verification; nonlinear moving thermal/energy and formation remain incomplete.',
        provenance=frozen,rounded_provenance=old,rounded_changed_sources=changed,summary=summary,
        cases=[compact(r) for r in new_records],rounded_cases=[compact(r) for r in old_records],audits=audits,rounded_audits=old_audits,
        previous_linear_diagnosis=previous,previous_diagnosis_sha256=digest(previous_path))
    files=helpers();snapshot=directory/'report_source_snapshot.tar.gz'
    with tarfile.open(snapshot,'w:gz') as tar:
        for name in files:tar.add(name,arcname=name)
    with tarfile.open(snapshot) as tar:
        for name,sha in files.items():assert hashlib.sha256(tar.extractfile(name).read()).hexdigest()==sha
    report.update(helper_sources=files,report_source_archive_sha256=digest(snapshot))
    figure=plot(report,directory);report['figures']={figure.name:digest(figure)}
    write_record(directory/'report.json',report);print(json.dumps(summary,indent=2),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/compatible_wave_verification')
    run(parser.parse_args().directory)
