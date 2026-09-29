"""Audit the complete smooth reference, preserve failed attempts and compare limits."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import tarfile
import numpy as np
from .smooth_finite_mass_cases import rate_protocol,evolution_protocol,initial
from .smooth_finite_mass_audit import source,rate,trajectory
from .finite_mass_audit import digest,load,state_from
from .finite_mass_time_comparison import measure as time_comparison
from .smooth_finite_mass_plots import plot
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
                if name.startswith(('examples.','plasma_bench.')):visit(name)
    for name in ('plasma_bench.__init__','examples.smooth_finite_mass_report','examples.smooth_precision_audit',
                 'examples.smooth_slow_mode_diagnosis','tests.test_smooth_finite_mass','tests.test_kernel_induction'):visit(name)
    paths=[Path(name.replace('.','/')+'.py') for name in sorted(seen)]
    paths+=[path.with_suffix('.md') for path in paths if path.with_suffix('.md').is_file()]
    return {str(path):digest(path) for path in sorted(paths)}


def initial_attempt(directory,current):
    record=json.loads((directory/'attempt.json').read_text());old=record['source']
    assert digest(directory/'source_snapshot.tar.gz')==record['source_archive_sha256']
    with tarfile.open(directory/'source_snapshot.tar.gz') as tar:
        for path,sha in old['source_files'].items():
            data=tar.extractfile(path).read();assert hashlib.sha256(data).hexdigest()==sha
            if path=='examples/smooth_finite_mass_measurements.py':
                expected=data.decode().replace('completed=failure is None and model.time==end,','completed=bool(failure is None and model.time==end),')
                assert Path(path).read_text()==expected
            else:assert digest(path)==sha,path
    matched=[]
    for kind in ('rates','trajectories'):
        for path in sorted((directory/kind).glob('*.json')):
            a,aa=load(directory/kind,path.stem);b,bb=load(current/kind,path.stem)
            assert set(aa)==set(bb)
            for key in aa:np.testing.assert_array_equal(aa[key],bb[key])
            for key in a:
                if key not in ('provenance','elapsed_s'):assert a[key]==b[key],(path.name,key)
            matched.append(str(path.relative_to(directory)))
    orphans={}
    for name in record['orphan_field_archives']:
        with np.load(directory/name) as a,np.load(current/name) as b:
            assert set(a.files)==set(b.files)
            for key in a.files:np.testing.assert_array_equal(a[key],b[key])
        orphans[name]='Exact fields reproduced; missing historical histories remain unverified.'
    return dict(manifest_sha256=digest(directory/'attempt.json'),source_sha256=old['source_sha256'],
        exact_reproduced_records=matched,orphan_field_replays=orphans)


def run(directory):
    directory=Path(directory);frozen=source(directory);rates=[];rate_audits={};runs=[];run_audits={}
    for case in rate_protocol():
        r,a=rate(directory,case,frozen['provenance']);rates.append(r);rate_audits[case['name']]=a
        print('rate audit',case['name'],r['completed'],flush=True)
    for case in evolution_protocol():
        r,a=trajectory(directory,case,frozen['provenance']);runs.append(r);run_audits[case['name']]=a
        print('trajectory audit',case['name'],r['completed'],flush=True)
    previous=json.loads(Path('outputs/finite_mass_evolution/report.json').read_text())
    for row in previous['rate_cases']:
        actual,arrays=load(Path('outputs/finite_mass_initial_attempt'),row['configuration']['name']);assert actual['metrics']==row['metrics']
    for row in previous['evolution_cases']:
        actual,arrays=load(Path('outputs/finite_mass_evolution'),row['configuration']['name']);assert actual['metrics']==row['metrics']
    earlier=initial_attempt(Path('outputs/smooth_finite_mass_initial_attempt'),directory)
    compact=lambda row:{key:value for key,value in row.items() if key not in ('history_columns','provenance','elapsed_s')}
    precision=json.loads((directory/'precision/contact_kernel.json').read_text())
    assert precision['current_source']==frozen['provenance']
    assert digest(directory/'precision'/precision['fields']['file'])==precision['fields']['sha256']
    slow=json.loads((directory/'slow_mode_diagnosis.json').read_text())
    assert digest(directory/slow['fields']['file'])==slow['fields']['sha256']
    summary=dict(rate_cases=len(rates),completed_rates=sum(bool(r['completed']) for r in rates),
        trajectories=len(runs),completed_trajectories=sum(bool(r['completed']) for r in runs),
        maximum_energy_budget=max(a.get('max_relative_energy_budget',0.) for a in run_audits.values()),
        maximum_integral_drift=max((value for a in run_audits.values() for value in a.get('normalized_integral_drifts',{}).values()),default=0.),
        maximum_norm_quadrature_change=max(a.get('norm_quadrature_absolute_change',0.) for a in [*rate_audits.values(),*run_audits.values()]),
        failures={r['configuration']['name']:r['failure'] for r in [*rates,*runs] if not r['completed']})
    report=dict(scope='Smooth finite-mass transport improvement; full magnetic, thermal, vacuum and formation accuracy remains incomplete.',
        provenance=frozen['provenance'],summary=summary,rates=[compact(r) for r in rates],rate_audits=rate_audits,
        trajectories=[compact(r) for r in runs],trajectory_audits=run_audits,
        baseline_report_sha256=digest('outputs/finite_mass_evolution/report.json'),baseline_rates=previous['rate_cases'],
        baseline_trajectories=previous['evolution_cases'],initial_attempt=earlier,
        precision=precision,slow_mode_diagnosis=slow,time_comparisons=time_comparison(directory/'trajectories'))
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
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/smooth_finite_mass_verification')
    run(parser.parse_args().directory)
