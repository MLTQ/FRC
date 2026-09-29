"""Audit all moving-current records and report independent physical evidence."""
from pathlib import Path
import json
import numpy as np
from .moving_current_measurements import protocol
from .moving_current_records import provenance
from .moving_current_audit import audit_case,audit_initial,digest
from .moving_species_audit import audit as audit_species
from .moving_current_plot import plot
from .verification_records import write_record


def run(directory=Path('outputs/moving_current_verification')):
    source=provenance();names=('examples/verify_moving_current.py','examples/moving_current_audit.py',
        'examples/moving_species_audit.py','examples/moving_current_plot.py','tests/test_moving_current.py','tests/test_moving_species.py')
    helpers={n:digest(Path(n)) for n in names};records=[];audits={}
    for case in protocol():
        record,result=audit_case(directory,case,source);records.append(record);audits[case['name']]=result
        print(case['name'],record['status'],'archive/physics audit passed',flush=True)
    lookup={r['configuration']['name']:r for r in records};refinement={}
    for prefix,grids in (('moving',(8,16,32,64,128)),('target',(8,16,32,64)),('massless',(8,16,32,64))):
        refinement[prefix]={}
        group=[lookup[f'{prefix}_{nr}'] for nr in grids]
        assert all(r['completed'] for r in group)
        for key in ('magnetic_rate','current_rate','electric'):
            refinement[prefix][key]={}
            for metric in ('relative_L2','axis_relative_L2'):
                errors=[r['metrics'][key][metric] for r in group]
                assert all(b<a for a,b in zip(errors,errors[1:])),(prefix,key,metric,errors)
                refinement[prefix][key][metric]=dict(radial_cells=list(grids),errors=errors,
                    observed_orders=[float(np.log2(a/b)) for a,b in zip(errors,errors[1:])])
    controls={}
    for control,base in (('quadrature_16','moving_16'),('tolerance_16','moving_16'),('tolerance_64','moving_64')):
        a,b=lookup[control],lookup[base];assert a['completed'] and b['completed']
        controls[control]={key:{metric:abs(a['metrics'][key][metric]-b['metrics'][key][metric])
                              for metric in ('relative_L2','axis_relative_L2')} for key in a['metrics']}
    species=audit_species();initial=audit_initial(Path('outputs/moving_current_initial_attempt'))
    plot(records,directory/'refinement.png')
    successful=[r for r in records if r['completed']]
    summary=dict(expected_cases=len(protocol()),completed=len(successful),rejected=len(records)-len(successful),
        max_linear_residual=max(r['linear_residual'] for r in successful),
        cases_using_linear_refinement=[r['configuration']['name'] for r in successful if r['linear_refinements']],
        max_full_electric_faraday_relative=max(a['full_electric_faraday_relative'] for a in audits.values() if 'full_electric_faraday_relative' in a),
        max_physical_energy_derivative_residual=max(d['relative_balance'] for a in audits.values() if 'energy_work' in a for d in a['energy_work']['independent_physical_energy_derivative']),
        max_norm_quadrature_absolute_change=max(v for a in audits.values() if 'norm_quadrature_absolute_change' in a for row in a['norm_quadrature_absolute_change'].values() for v in row.values()))
    assert source==provenance() and helpers=={n:digest(Path(n)) for n in names}
    report=dict(audit_passed=True,provenance=source,helper_sources=helpers,summary=summary,cases=records,audits=audits,
        refinement=refinement,quadrature_and_tolerance_controls=controls,parent_species=species,initial_incomplete_attempt=initial,
        figure=dict(file='refinement.png',sha256=digest(directory/'refinement.png')),
        scope='Instantaneous positive-density bounded manufactured induction rates. Parent-species continuum stress/energy contract; no conservative moving time integration, physical free edge or formation validation.')
    write_record(directory/'report.json',report);print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':run()
