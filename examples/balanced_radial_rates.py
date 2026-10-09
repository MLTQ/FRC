"""Full physical rates at fixed independent radial wavelength on separate grids."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.balanced_wall_plasma import balanced
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .matched_mode_rates import measure
from .radial_mode_interfaces import measure as interfaces
from .implicit_plasma_records import provenance,freeze,sha


def calculate(spec,mode,orders=(16,20)):
    start=time.monotonic();epsilon=spec.get('epsilon',1e-4)
    if not np.isfinite(epsilon) or epsilon<=0:raise ValueError('Positive finite perturbation amplitude required')
    d,u,system,bias,du,db=fixture(spec,mode);system=balanced(system)
    arrays=dict(state=u,fluid_perturbation=du,magnetic_perturbation=db);rates=[]
    result=dict(spec=spec,epsilon=epsilon,input_sha256=mode.input_sha256,input_summary_sha256=mode.input_summary_sha256,
        radial_mode=mode.radial_mode,reference_points=mode.points,continuum_growth=mode.growth)
    for name,sign in [('minus',-1),('plus',1)]:
        state=u+sign*epsilon*du;field=field_from_packed(system.geometry,sign*epsilon*db,bias.axial_bias)
        arrays.update({name+'_state':state,name+'_magnetic':field.packed_increment(system.geometry)})
        try:rate=system.derivative(state,field)
        except (RuntimeError,ValueError) as error:
            result.update(status='failed',failed_sign=name,failure_type=type(error).__name__,failure_message=str(error),elapsed_seconds=time.monotonic()-start)
            return result,arrays
        rates.append(rate)
        arrays.update({name+'_'+key:value for key,value in dict(fluid_rate=rate.fluid,magnetic_rate=rate.magnetic,
            electric=rate.electric,current_rate=rate.current).items()})
    quadrature=[]
    for order in orders:
        physical=measure(mode,d,system,u,du,db,epsilon,rates,order)
        face,raw=interfaces(mode,d,system,u,db,order);physical['representation_and_faces']=face
        arrays.update({f'q{order}_{key}':value for key,value in raw.items()});quadrature.append(physical)
    result.update(status='complete',quadrature=quadrature,max_induction_residual=max(r.residual for r in rates),
        rate_diagnostics=[r.diagnostics for r in rates],elapsed_seconds=time.monotonic()-start)
    return result,arrays


def run(destination,specs,reference_directory='outputs/radial16_reference',reference_index=1):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    mode=ArchivedContactMode(reference_directory,reference_index)
    identity=provenance(['examples.balanced_radial_rates']);destination.mkdir(parents=True);freeze(destination,identity)
    records=[];summary=dict(terminal=False,source_sha256=identity['source_sha256'],cases=records,
        reference_directory=str(reference_directory),reference_index=reference_index)
    for spec in specs:
        record,arrays=calculate(spec,mode);path=destination/f'case_{len(records):02}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
        (destination/f'case_{len(records)-1:02}.json').write_text(json.dumps(record,indent=2)+'\n')
        (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        metrics=record.get('quadrature',[{}])[-1].get('rate_errors',{})
        print('fixed radial mode',spec,record['status'],{k:metrics[k]['relative_L2'] for k in ('radial_velocity','radial_B','current_rate') if k in metrics},flush=True)
    if provenance(['examples.balanced_radial_rates'])!=identity:raise RuntimeError('Source changed')
    summary.update(terminal=True,all_cases_complete=all(r['status']=='complete' for r in records))
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
