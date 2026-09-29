"""Archive and independently audit current-increment precision controls."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from .current_thermal_cases import protocol
from .current_thermal_evolution import measure
from .current_thermal_records import provenance,save
from .current_thermal_audit import audit_case
from .current_stationary_measurements import measure as stationary
from .verification_records import write_record


def cases():
    base={c['name']:c for c in protocol()};result=[]
    for nr in (8,16,32):
        result.extend(deepcopy(base[name]) for name in (f'original_{nr}',f'wave_target_{nr}'))
        for amplitude in (1e-9,1e-12):
            case=deepcopy(base[f'wave_target_{nr}'])
            case['name']=f'weak_{amplitude:g}_{nr}';case['reference']['amplitude']=amplitude
            result.append(case)
    for label,changes in (('strong_bias',dict(bias=.5)),('reverse',dict(helicity=-1)),
                           ('strong_inertia',dict(number_density=1e18))):
        case=deepcopy(base['wave_target_16']);case['name']='weak_'+label
        case['reference'].update(amplitude=1e-12,**changes);result.append(case)
    return result


def run(directory=Path('outputs/current_precision_verification')):
    directory.mkdir(parents=True,exist_ok=True);source=provenance()
    names=('examples/verify_current_precision.py','examples/current_stationary_measurements.py',
           'examples/current_thermal_audit.py','tests/test_current_precision.py')
    helpers={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}
    stationary_cases=[];records=[];audits={}
    for nr in (8,16,32):
        for bias in (.005,-.005,.5):
            for dt in (1e-8,1e-6):
                record=stationary(directory,nr,bias,dt);stationary_cases.append(record)
                print(record['name'],max(h['magnetic_max_error_t'] for h in record['history']),flush=True)
    for case in cases():
        record_path=directory/(case['name']+'.json')
        if record_path.exists():
            record=json.loads(record_path.read_text())
            assert record['configuration']==case and record['provenance']==source
            assert hashlib.sha256((directory/record['fields']['file']).read_bytes()).hexdigest()==record['fields']['sha256']
        else:
            record,arrays=measure(case);record=save(directory,case,source,record,arrays)
        assert record['completed'],(case['name'],record['failure'])
        with np.load(directory/record['fields']['file'],allow_pickle=False) as z:
            archived={k:z[k] for k in z.files}
        mesh=MHDMesh(case['reference']['radius'],case['reference']['length']/2,case['nr'],case['nz'])
        audits[case['name']]=audit_case(case,record,archived,mesh);records.append(record)
        print(case['name'],record['metrics'],flush=True)
    lookup={r['configuration']['name']:r for r in records};refinement={}
    for amplitude in (1e-9,1e-12):
        group=[lookup[f'weak_{amplitude:g}_{nr}'] for nr in (8,16,32)]
        for key in ('magnetic_relative_L2','current_relative_L2'):
            errors=[r['metrics'][key] for r in group]
            assert errors[1]<.6*errors[0] and errors[2]<.6*errors[1]
            refinement[f'{amplitude:g}_{key}']=errors
    comparisons={}
    for record in records:
        name=record['configuration']['name'];old=Path('outputs/current_thermal_verification')/(name+'.json')
        if old.exists():
            baseline=json.loads(old.read_text())
            comparisons[name]=dict(baseline_source=baseline['provenance']['source_sha256'],
                baseline_file_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
                metric_absolute_changes={k:record['metrics'][k]-v for k,v in baseline['metrics'].items()})
    assert source==provenance() and helpers=={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}
    summary=dict(stationary_cases=len(stationary_cases),evolution_cases=len(records),
        max_stationary_relative_B_error=max(h['magnetic_max_error_t']/abs(c['bias_t']) for c in stationary_cases for h in c['history']),
        max_energy_relative_drift=max(a['max_total_energy_j_relative_drift'] for a in audits.values()),
        max_stored_state_residual=max(max(r['history_columns']['stored_state_residual']) for r in records),
        max_nonlinear_residual=max(max(r['history_columns']['nonlinear_residual']) for r in records))
    report=dict(audit_passed=True,provenance=source,helper_sources=helpers,summary=summary,
        stationary=stationary_cases,cases=records,audits=audits,weak_refinement=refinement,baseline_comparisons=comparisons,
        scope='Fixed positive density and constrained zero barycentric flow; precision fix, not moving or free-edge formation validation.')
    write_record(directory/'report.json',report);print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':run()
