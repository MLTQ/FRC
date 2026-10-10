"""Frozen controls for shared driving coefficients and term quadrature."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    first=Path('outputs/radial16_term_capture/case_01.npz');second=Path('outputs/radial16_drive_capture/case_00.npz')
    a=dict(np.load(first));b=dict(np.load(second));keys=sorted(a.keys()&b.keys())
    if not keys or any(not np.array_equal(a[k],b[k]) for k in keys):raise ValueError('Shared baseline coefficients changed')
    records=[];inputs={str(p):sha(p) for p in (first,second)}
    for directory in ('radial16_term_audit','radial16_drive_audit'):
        path=Path('outputs')/directory/'summary.json';summary=json.loads(path.read_text())
        if not summary['terminal']:raise ValueError('Terminal audits required')
        inputs[str(path)]=sha(path)
        for i,case in enumerate(summary['cases']):
            low,high=case['quadrature'];differences=[]
            for term,regions in low['metrics'].items():
                for region,values in regions.items():
                    for name in ('actual_L2','target_L2','error_L2'):
                        x=np.array(values[name]);y=np.array(high['metrics'][term][region][name]);scale=np.maximum(np.maximum(abs(x),abs(y)),np.finfo(float).tiny)
                        for component in range(3):
                            differences.append(dict(term=term,region=region,norm=name,component=component,relative=float(abs(x[component]-y[component])/scale[component]),absolute=float(abs(x[component]-y[component])),low=float(x[component]),high=float(y[component])))
            physical=[v for v in differences if v['region']=='volume']
            records.append(dict(directory=directory,index=i,largest_relative=max(differences,key=lambda x:x['relative']),
                largest_volume_relative=max(physical,key=lambda x:x['relative']),largest_absolute=max(differences,key=lambda x:x['absolute'])))
    path=Path('outputs/radial16_roundoff_probe/summary.json');probe=json.loads(path.read_text());inputs[str(path)]=sha(path)
    if not probe['terminal'] or not probe['matches_original_first_attempt_failure']:raise ValueError('Exact original rejection required')
    identity=provenance(['examples.radial_diagnostic_comparison']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_hashes=inputs,
        baseline_shared_coefficients_bitwise_equal=True,shared_coefficient_keys=keys,quadrature_controls=records,
        original_first_attempt_failure_bitwise_reproduced=True,failed_full_case_remains_failed=True)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
