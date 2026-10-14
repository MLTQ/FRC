"""Unchanged complete spatial RHS after verified fluid-direction subdivision."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.regular_current import RegularCurrent
from .split_mode_audit import fixture
from .trajectory_tangent_probe import verified_inputs,constraints
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,transfer_directory,elements):
    destination=Path(destination);transfer_directory=Path(transfer_directory);elements=int(elements)
    if destination.exists():raise ValueError('Fresh full refinement-probe output required')
    if elements not in (16,32):raise ValueError('Declared fine fluid grid required')
    transferred,inputs=verified_inputs(transfer_directory,str(transfer_directory)+'_replay_verification')
    inputs.update(transferred['input_sha256'])
    if not transferred['all_cases_complete'] or transferred['maximum_differential_relative']>transferred['transfer_relative_bound']:
        raise ValueError('Complete passing physical direction transfer required')
    grid=next(g for g in transferred['grids'] if g['elements']==elements)
    context_path=transfer_directory/grid['context_file']
    if sha(context_path)!=grid['context_sha256']:raise ValueError('Original transferred context required')
    with np.load(context_path) as raw:context=dict(raw)
    inputs[str(context_path)]=sha(context_path)
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original replayed inputs required')
    spec=dict(transferred['original_spec'],elements=elements)
    identity=provenance(['examples.coupled_fluid_probe']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,
                transfer_directory=str(transfer_directory),spec=spec,amplitudes=transferred['amplitudes'],cases=[],directions={})
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save()
    try:
        mode,d,base,system,bias,_,_=fixture(spec)
        if not np.array_equal(base,context['state']):raise ValueError('Original fixture base must reproduce exact transferred base')
        if bias.axial_bias!=transferred['physical_bias'] or system.solver_options['tolerance']!=1e-10:
            raise ValueError('Original bias and induction bound required')
        coordinates=ModalWallCoordinates(system,base,bias);model=RegularCurrent(system.mapping,base[...,0],system.ion_mass)
        context.update(fluid_row_scale=coordinates.fluid_row_scale,radial_basis=coordinates.radial_basis,
                       magnetic_scalars=coordinates.magnetic_scalars)
        for name in ('saved_error','independent_mode'):
            du,db=context[name+'_fluid'],context[name+'_magnetic']
            actual=model.current(db)
            if not np.array_equal(actual,context[name+'_current']):raise ValueError('Original current must remain bitwise identical')
            value=coordinates.pack(du,db);context[name+'_coordinate']=value
            result['directions'][name]=dict(normalization=transferred['directions'][name]['normalization'],
                coordinate_max=float(np.max(abs(value))),constraints=constraints(system,du,db),
                current_coefficients_exact=True,normalization_changed=False)
        path=destination/'context.npz';np.savez_compressed(path,**context)
        result.update(context_file=path.name,context_sha256=sha(path),continuum_growth=mode.growth,
            reference_input_sha256=mode.input_sha256,physical_bias=bias.axial_bias,
            induction_tolerance=system.solver_options['tolerance'],coordinate_size=coordinates.size);save()
        before=(fingerprint(base),fingerprint(bias.packed_increment(system.geometry)))
        base_rate=system.derivative(base,bias)
        if before!=(fingerprint(base),fingerprint(bias.packed_increment(system.geometry))):raise ValueError('Base input modified')
        if base_rate.residual>result['induction_tolerance']:raise ValueError('Base original induction bound failed')
        path=destination/'base_rate.npz'
        np.savez_compressed(path,fluid_rate=base_rate.fluid,magnetic_rate=base_rate.magnetic,electric=base_rate.electric,current_rate=base_rate.current)
        result['base']=dict(arrays_file=path.name,arrays_sha256=sha(path),residual=base_rate.residual,
            diagnostics=base_rate.diagnostics,actual_inputs_byte_unchanged=True,
            coordinate_rate_max=float(np.max(abs(coordinates.pack(base_rate.fluid,base_rate.magnetic)))));save()
        for name in ('saved_error','independent_mode'):
            du,db=context[name+'_fluid'],context[name+'_magnetic'];value=context[name+'_coordinate']
            for amplitude in result['amplitudes']:
                index=len(result['cases']);stem=f'case_{index:02}';started=time.monotonic()
                arrays={sign+'_'+key:a for sign,multiplier in [('minus',-1),('plus',1)]
                    for key,a in dict(state=base+multiplier*amplitude*du,magnetic=context['magnetic']+multiplier*amplitude*db).items()}
                path=destination/(stem+'_inputs.npz');np.savez_compressed(path,**arrays)
                row=dict(direction=name,amplitude=amplitude,status='running',input_file=path.name,input_sha256=sha(path),signs=[])
                result['cases'].append(row);save();rates=[]
                for sign in ('minus','plus'):
                    state=arrays[sign+'_state'];field=coordinates.field(arrays[sign+'_magnetic'])
                    before=(fingerprint(state),fingerprint(field.packed_increment(system.geometry)))
                    sign_record=dict(sign=sign,status='running');row['signs'].append(sign_record);save()
                    try:rate=system.derivative(state,field)
                    except Exception as error:
                        sign_record.update(status='rhs_exception',exception_type=type(error).__name__,message=str(error));row['status']='rhs_exception'
                    finally:
                        after=(fingerprint(state),fingerprint(field.packed_increment(system.geometry)))
                        sign_record['actual_inputs_byte_unchanged']=before==after;save()
                    if before!=after:raise ValueError('Actual passed state or field modified')
                    if row['status']=='rhs_exception':break
                    if rate.residual>result['induction_tolerance']:raise ValueError('Original induction bound failed')
                    path=destination/(stem+'_'+sign+'.npz')
                    np.savez_compressed(path,fluid_rate=rate.fluid,magnetic_rate=rate.magnetic,electric=rate.electric,current_rate=rate.current)
                    rates.append(rate);fields=d.fields(state)
                    sign_record.update(status='complete',arrays_file=path.name,arrays_sha256=sha(path),residual=rate.residual,
                        diagnostics=rate.diagnostics,point_density_min=float(state[...,0].min()),
                        point_pressure_min_by_species=np.min(fields['pressure'],axis=(0,1,2)).tolist(),
                        constraints=constraints(system,state,field.packed_increment(system.geometry)));save()
                if len(rates)==2:
                    fluid=(rates[1].fluid-rates[0].fluid)/(2*amplitude);magnetic=(rates[1].magnetic-rates[0].magnetic)/(2*amplitude)
                    action=coordinates.pack(fluid,magnetic);projection=float(np.dot(value,action)/np.dot(value,value))
                    path=destination/(stem+'_action.npz');np.savez_compressed(path,fluid=fluid,magnetic=magnetic,coordinate=action,current=model.current(magnetic))
                    action_norm=float(np.linalg.norm(action))
                    row.update(status='complete',action_file=path.name,action_sha256=sha(path),coordinate_projection=projection,
                        coordinate_orthogonal_relative=float(np.linalg.norm(action-projection*value)/action_norm) if action_norm else None)
                row['elapsed_seconds']=time.monotonic()-started;save()
                print('full coupled fluid refinement',elements,name,amplitude,row['status'],row.get('coordinate_projection'),flush=True)
        result['all_cases_complete']=all(c['status']=='complete' for c in result['cases'])
    except BaseException as error:
        result.update(exception_type=type(error).__name__,exception_message=str(error));raise
    finally:
        result['terminal']=True
        if provenance(['examples.coupled_fluid_probe'])!=identity or any(sha(p)!=h for p,h in inputs.items()):
            result['source_or_input_changed']=True;save();raise ValueError('Probe source/input changed')
        save()
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('transfer_directory');p.add_argument('elements')
    a=p.parse_args();run(a.destination,a.transfer_directory,a.elements)
