"""Audit and summarize the complete conservative finite-mass baseline study."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import tarfile
import numpy as np
from .finite_mass_audit import baseline,evolution,digest,load
from .finite_mass_evolution import protocol
from .finite_mass_records import provenance
from .finite_mass_flux_diagnosis import measure as flux_diagnosis
from .finite_mass_time_comparison import measure as time_comparison
from .finite_mass_plots import plot
from .verification_records import write_record


def helper_sources():
    """Archive the full local import closure, companions and focused tests."""
    seen=set()
    def visit(module):
        path=Path(module.replace('.','/')+'.py')
        if module in seen or not path.is_file():return
        seen.add(module)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.ImportFrom) and node.module:
                names=['.'.join(module.split('.')[:-node.level]+node.module.split('.')) if node.level else node.module]
            elif isinstance(node,ast.Import):names=[a.name for a in node.names]
            else:continue
            for name in names:
                if name.startswith(('plasma_bench.','examples.')):visit(name)
    for module in ('plasma_bench.__init__','examples.finite_mass_report','tests.test_finite_mass_mhd',
                   'tests.test_finite_mass_wave','tests.test_annular_polynomial'):visit(module)
    paths=[Path(module.replace('.','/')+'.py') for module in sorted(seen)]
    paths+=[path.with_suffix('.md') for path in paths if path.with_suffix('.md').is_file()]
    return {str(path):digest(path) for path in sorted(paths)}


def run(baseline_directory,directory):
    baseline_directory=Path(baseline_directory);directory=Path(directory);source=provenance()
    manifest,rates,rate_audits=baseline(baseline_directory)
    runs=[];audits={}
    for case in protocol():
        record,audit=evolution(directory,case,source);runs.append(record);audits[case['name']]=audit
        print('trajectory audit',case['name'],record['completed'],flush=True)
    diagnosis,arrays=flux_diagnosis();path=directory/'flux_diagnosis.npz';np.savez_compressed(path,**arrays)
    diagnosis['fields']=dict(file=path.name,sha256=digest(path))
    time=time_comparison(directory)
    compact=lambda row:{key:row[key] for key in ('configuration','parameters','completed','failure','metrics','fields')}
    mass_drift=max(max(abs(np.array(row['history_columns']['mass_kg'])/row['history_columns']['mass_kg'][0]-1)) for row in runs)
    unforced=[]
    for row in runs:
        if row['configuration']['family']!='unforced':continue
        h=row['history_columns'];name=row['configuration']['name']
        unforced.append(dict(name=name,electron_physical_entropy_change_per_particle=audits[name]['electron_physical_entropy_change_per_particle'],
            ion_physical_entropy_change_per_particle=audits[name]['ion_physical_entropy_change_per_particle']))
    report=dict(scope='Conservative positive-density moving finite-mass reference; temperature accuracy is not adequate for formation. No vacuum or transport closure validation.',
        provenance=source,baseline_provenance=manifest['provenance'],baseline_manifest_sha256=digest(baseline_directory/'baseline.json'),
        rate_cases=[compact(row) for row in rates],rate_audits=rate_audits,
        evolution_cases=[compact(row) for row in runs],evolution_audits=audits,
        flux_diagnosis=diagnosis,time_comparisons=time,
        summary=dict(rate_cases=len(rates),completed_rates=sum(row['completed'] for row in rates),
            trajectories=len(runs),completed_trajectories=sum(row['completed'] for row in runs),
            maximum_relative_mass_drift=float(mass_drift),maximum_relative_energy_budget=max(a['max_relative_energy_budget'] for a in audits.values()),
            minimum_electron_pressure_pa=min(min(row['history_columns']['min_electron_pressure_pa']) for row in runs),
            minimum_ion_pressure_pa=min(min(row['history_columns']['min_ion_pressure_pa']) for row in runs),
            maximum_norm_quadrature_change=max([a['norm_quadrature_absolute_change'] for a in audits.values()]+[a['norm_quadrature_absolute_change'] for a in rate_audits.values()]),
            unforced_entropy=unforced),
        limitations=['Slow nearly solenoidal mass/pressure rates remain very inaccurate under the limited central reconstruction.',
            'Fast-wave and manufactured ion pressure remain much less accurate than conserved global energy.',
            'Cubic continuity reconstruction is diagnostic only; full consistent momentum/enthalpy/kinetic coupling is required.',
            'No entropy inequality, shock heating partition, stable general timestep or plasma-vacuum formation claim.'])
    assert source==provenance()
    helpers=helper_sources();archive=directory/'source_snapshot.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for name in helpers:tar.add(name,arcname=name)
    with tarfile.open(archive) as tar:
        for name,sha in helpers.items():assert hashlib.sha256(tar.extractfile(name).read()).hexdigest()==sha
    report.update(helper_source_files=helpers,source_archive_sha256=digest(archive))
    figures=plot(report,directory);report['figures']={path.name:digest(path) for path in figures}
    write_record(directory/'report.json',report)
    print(json.dumps(report['summary'],indent=2),flush=True)
    print('time comparisons',json.dumps(time,indent=2),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--baseline-directory',default='outputs/finite_mass_initial_attempt')
    parser.add_argument('--directory',default='outputs/finite_mass_evolution');args=parser.parse_args()
    run(args.baseline_directory,args.directory)
