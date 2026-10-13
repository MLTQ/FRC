"""Strict full spatial derivatives along a saved error and an independent mode."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.regular_current import RegularCurrent
from .split_mode_audit import fixture
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance, freeze, sha


def verified_inputs(directory, verification):
    directory=Path(directory);verification=Path(verification)
    original=json.loads((directory/'summary.json').read_text())
    replay=json.loads((verification/'summary.json').read_text())
    source=json.loads((directory/'provenance.json').read_text())
    if not original['terminal'] or not replay['terminal'] or not replay['bitwise_metrics']:
        raise ValueError('Terminal trajectory and complete protocol replay required')
    if replay['original_directory']!=str(directory):raise ValueError('Matching original trajectory required')
    if any(sha(p)!=h for p,h in source['source_files'].items()):raise ValueError('Original physical source required')
    inputs={str(directory/'provenance.json'):sha(directory/'provenance.json'),
        str(verification/'summary.json'):sha(verification/'summary.json')}
    for prefix in ('original','repeated'):
        p=Path(replay[prefix+'_directory'])/'summary.json'
        if sha(p)!=replay[prefix+'_summary_sha256']:raise ValueError('Unchanged complete replay required')
        inputs[str(p)]=sha(p)
    for name,row in replay['arrays'].items():
        if row['original_sha256']!=row['repeated_sha256']:raise ValueError('Byte-identical full protocol required')
        for prefix in ('original','repeated'):
            p=Path(replay[prefix+'_directory'])/name
            if sha(p)!=row[prefix+'_sha256']:raise ValueError('Original replayed raw files required')
            inputs[str(p)]=sha(p)
    return original,inputs


def constraints(system,state,magnetic):
    psi,g=system.geometry.unpack(magnetic)
    return dict(axis_transverse_max=float(np.max(abs(state[0,0,:,1:3]))),
        wall_normal_momentum_max=float(np.max(abs(state[-1,-1,:,1]))),
        fixed_psi_max=float(np.max(abs(psi[[0,-1]]))),
        toroidal_wall_trace_max=float(np.max(abs(system.geometry.wall_trace@g))))


def run(destination, trajectory, verification):
    destination=Path(destination);trajectory=Path(trajectory)
    if destination.exists():raise ValueError('Fresh tangent-probe destination required')
    original,input_hashes=verified_inputs(trajectory,verification)
    if original['spec']!={'elements':8,'points':64,'field_nr':32,'field_nz':64} or original['completed_steps']<10:
        raise ValueError('Declared refined contact with accepted step10 required')
    seeds=['examples.trajectory_tangent_probe','examples.split_mode_audit',
        'examples.matched_mode_fixture','examples.matched_mode_audit','examples.matched_mode_comparison']
    identity=provenance(seeds);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],
        trajectory=str(trajectory),verification=str(verification),input_sha256=input_hashes,
        spec=original['spec'],selected_step=10,amplitudes=[1e-4,1e-5,1e-6],cases=[])
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save()
    try:
        mode,d,base,system,bias,reference_u,reference_b=fixture(original['spec'])
        if mode.input_sha256!=original['input_sha256'] or mode.growth!=original['continuum_growth']:
            raise ValueError('Original independent continuum input required')
        tolerance=system.solver_options['tolerance']
        if tolerance!=1e-10:raise ValueError('Original induction bound required')
        coordinates=ModalWallCoordinates(system,base,bias)
        with np.load(trajectory/'initial.npz') as raw:initial={k:raw[k].copy() for k in raw.files}
        selected=original['records'][9]
        with np.load(trajectory/selected['arrays_file']) as raw:accepted={k:raw[k].copy() for k in raw.files}
        expected_u=base+original['epsilon']*reference_u
        expected_b=coordinates.field(original['epsilon']*reference_b).packed_increment(system.geometry)
        if fingerprint(initial['state'])!=fingerprint(expected_u) or fingerprint(initial['magnetic'])!=fingerprint(expected_b):
            raise ValueError('Exact original initialization required')
        factor=original['epsilon']*np.expm1(mode.growth*selected['time'])
        error_u=(accepted['state']-initial['state'])-factor*reference_u
        error_b=(accepted['magnetic']-initial['magnetic'])-factor*reference_b
        zero_b=bias.packed_increment(system.geometry)
        if np.any(zero_b):raise ValueError('Original stationary zero magnetic increment required')
        model=RegularCurrent(system.mapping,base[...,0],system.ion_mass)
        context=dict(state=base,magnetic=zero_b,reference_fluid=reference_u,reference_magnetic=reference_b,
            native_error_fluid=error_u,native_error_magnetic=error_b,
            fluid_row_scale=coordinates.fluid_row_scale,radial_basis=coordinates.radial_basis,
            magnetic_scalars=coordinates.magnetic_scalars)
        directions={};metadata={}
        for name,fluid,magnetic in [('saved_error',error_u,error_b),('independent_mode',reference_u,reference_b)]:
            packed=coordinates.pack(fluid,magnetic);normalization=float(np.max(abs(packed)))
            if not np.isfinite(normalization) or normalization<=0:raise ValueError('Finite nonzero complete tangent required')
            value=packed/normalization;du,db=coordinates.unpack(value);directions[name]=(du,db,value)
            count=(system.geometry.mesh.nr+1)*system.geometry.mesh.nz
            roundtrip_u=du*normalization-fluid;roundtrip_b=db*normalization-magnetic
            metadata[name]=dict(normalization=normalization,
                pack_unpack_coordinate_max=float(np.max(abs(coordinates.pack(du,db)-value))),
                native_fluid_roundtrip_max_by_component=np.max(abs(roundtrip_u),axis=(0,1,2)).tolist(),
                native_magnetic_roundtrip_max_psi=float(np.max(abs(roundtrip_b[:count]))),
                native_magnetic_roundtrip_max_g=float(np.max(abs(roundtrip_b[count:]))),
                constraints=constraints(system,du,db))
            context.update({name+'_'+k:v for k,v in dict(fluid=du,magnetic=db,coordinate=value,current=model.current(db)).items()})
        path=destination/'context.npz';np.savez_compressed(path,**context)
        result.update(context_file=path.name,context_sha256=sha(path),directions=metadata,
            continuum_growth=mode.growth,reference_input_sha256=mode.input_sha256,
            original_epsilon=original['epsilon'],selected_time=selected['time'],coordinate_size=coordinates.size,
            induction_tolerance=tolerance,physical_bias=bias.axial_bias);save()
        base_before=(fingerprint(base),fingerprint(bias.packed_increment(system.geometry)))
        base_rate=system.derivative(base,bias)
        base_after=(fingerprint(base),fingerprint(bias.packed_increment(system.geometry)))
        if base_before!=base_after:raise ValueError('Actual base input modified')
        base_arrays={k:getattr(base_rate,v) for k,v in [('fluid_rate','fluid'),('magnetic_rate','magnetic'),('electric','electric'),('current_rate','current')]}
        path=destination/'base_rate.npz';np.savez_compressed(path,**base_arrays)
        result['base']=dict(arrays_file=path.name,arrays_sha256=sha(path),residual=base_rate.residual,
            diagnostics=base_rate.diagnostics,actual_inputs_byte_unchanged=True,
            coordinate_rate_max=float(np.max(abs(coordinates.pack(base_rate.fluid,base_rate.magnetic)))))
        if base_rate.residual>tolerance:raise ValueError('Base original bound failed')
        save()
        for name,(du,db,value) in directions.items():
            for amplitude in result['amplitudes']:
                started=time.monotonic();index=len(result['cases']);stem=f'case_{index:02}'
                arrays={sign+'_'+key:a for sign,multiplier in [('minus',-1),('plus',1)]
                    for key,a in dict(state=base+multiplier*amplitude*du,magnetic=zero_b+multiplier*amplitude*db).items()}
                path=destination/(stem+'_inputs.npz');np.savez_compressed(path,**arrays)
                row=dict(direction=name,amplitude=amplitude,status='running',input_file=path.name,input_sha256=sha(path),signs=[])
                result['cases'].append(row);save();rates=[]
                for sign in ('minus','plus'):
                    state=arrays[sign+'_state'];field=coordinates.field(arrays[sign+'_magnetic'])
                    before=(fingerprint(state),fingerprint(field.packed_increment(system.geometry)))
                    sign_record=dict(sign=sign,status='running');row['signs'].append(sign_record);save()
                    try:rate=system.derivative(state,field)
                    except Exception as error:
                        sign_record.update(status='rhs_exception',exception_type=type(error).__name__,message=str(error))
                        row['status']='rhs_exception'
                    finally:
                        after=(fingerprint(state),fingerprint(field.packed_increment(system.geometry)))
                        sign_record['actual_inputs_byte_unchanged']=before==after
                        save()
                    if before!=after:raise ValueError('Actual passed state or field modified')
                    if row['status']=='rhs_exception':break
                    if rate.residual>tolerance:raise ValueError('Original induction bound failed')
                    path=destination/(stem+'_'+sign+'.npz')
                    np.savez_compressed(path,fluid_rate=rate.fluid,magnetic_rate=rate.magnetic,electric=rate.electric,current_rate=rate.current)
                    rates.append(rate);fields=d.fields(state)
                    sign_record.update(status='complete',arrays_file=path.name,arrays_sha256=sha(path),residual=rate.residual,
                        diagnostics=rate.diagnostics,point_density_min=float(state[...,0].min()),
                        point_pressure_min_by_species=np.min(fields['pressure'],axis=(0,1,2)).tolist(),
                        constraints=constraints(system,state,field.packed_increment(system.geometry)));save()
                if len(rates)==2:
                    fluid=(rates[1].fluid-rates[0].fluid)/(2*amplitude)
                    magnetic=(rates[1].magnetic-rates[0].magnetic)/(2*amplitude)
                    action=coordinates.pack(fluid,magnetic)
                    projection=float(np.dot(value,action)/np.dot(value,value))
                    path=destination/(stem+'_action.npz')
                    np.savez_compressed(path,fluid=fluid,magnetic=magnetic,coordinate=action,current=model.current(magnetic))
                    action_norm=float(np.linalg.norm(action))
                    row.update(status='complete',action_file=path.name,action_sha256=sha(path),coordinate_projection=projection,
                        coordinate_orthogonal_relative=float(np.linalg.norm(action-projection*value)/action_norm) if action_norm else None)
                row['elapsed_seconds']=time.monotonic()-started;save()
                print('full spatial tangent',name,amplitude,row['status'],row.get('coordinate_projection'),flush=True)
        result['all_cases_complete']=all(c['status']=='complete' for c in result['cases'])
    except BaseException as error:
        result.update(exception_type=type(error).__name__,exception_message=str(error));raise
    finally:
        result['terminal']=True
        if provenance(seeds)!=identity or any(sha(p)!=h for p,h in input_hashes.items()):
            result['source_or_input_changed']=True;save();raise ValueError('Probe source or inputs changed')
        save()
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('trajectory');parser.add_argument('verification')
    args=parser.parse_args();run(args.destination,args.trajectory,args.verification)
