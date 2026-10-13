"""Independent physical attribution, continuum controls and interior comparisons."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler
from .trajectory_tangent_audit import expand, components, compare, NAMES
from .trajectory_density_prefix import fingerprint
from .trajectory_spatial_parts import GROUPS
from .spatial_parts_reference import analytic, interior
from .implicit_plasma_records import provenance, freeze, sha


def values(sample, states, rates, amplitude, prefix=''):
    minus=expand(sample,states['minus'],rates['minus'][prefix+'fluid'],rates['minus'][prefix+'magnetic'])
    plus=expand(sample,states['plus'],rates['plus'][prefix+'fluid'],rates['plus'][prefix+'magnetic'])
    fields=(plus[0]-minus[0])/(2*amplitude);pressure=(plus[1]-minus[1])/(2*amplitude)
    current=sample.electric((rates['plus'][prefix+'current']-rates['minus'][prefix+'current'])/(2*amplitude))
    electric=sample.electric((rates['plus'][prefix+'electric']-rates['minus'][prefix+'electric'])/(2*amplitude))
    return components(fields,pressure,current,electric),dict(rhs_fields=fingerprint(fields),
            rhs_pressure=fingerprint(pressure),current_rate=fingerprint(current),electric=fingerprint(electric))


def run(destination,parts_directory):
    destination=Path(destination);parts_directory=Path(parts_directory)
    if destination.exists():raise ValueError('Fresh physical contribution audit required')
    record=json.loads((parts_directory/'summary.json').read_text())
    source=json.loads((parts_directory/'provenance.json').read_text())
    if not record['terminal'] or not record['all_cases_complete'] or not record['all_recombined_original_bounds_pass']:
        raise ValueError('Complete original bounded contribution protocol required')
    if any(sha(p)!=h for p,h in source['source_files'].items()):raise ValueError('Original contribution source required')
    verifier=Path(str(parts_directory)+'_replay_verification')/'summary.json'
    replay=json.loads(verifier.read_text())
    if not replay['terminal'] or not replay['bitwise_metrics'] or replay['original_directory']!=str(parts_directory):
        raise ValueError('Complete matching contribution replay required')
    inputs={str(verifier):sha(verifier),str(parts_directory/'provenance.json'):sha(parts_directory/'provenance.json')}
    for side in ('original','repeated'):
        root=Path(replay[side+'_directory']);p=root/'summary.json'
        if sha(p)!=replay[side+'_summary_sha256']:raise ValueError('Unchanged complete replay summary required')
        inputs[str(p)]=sha(p)
        for name,row in replay['arrays'].items():
            if row['original_sha256']!=row['repeated_sha256'] or sha(root/name)!=row[side+'_sha256']:
                raise ValueError('Unchanged complete replay arrays required')
            inputs[str(root/name)]=sha(root/name)
    inputs.update(record['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original probe inputs required')
    probe=Path(record['probe']);probe_record=json.loads((probe/'summary.json').read_text())
    audit_path=Path('outputs/trajectory_tangent_audit/summary.json');old_audit=json.loads(audit_path.read_text())
    inputs[str(audit_path)]=sha(audit_path)
    audit_verifier=Path('outputs/trajectory_tangent_audit_replay_verification/summary.json')
    old_verify=json.loads(audit_verifier.read_text());inputs[str(audit_verifier)]=sha(audit_verifier)
    if not old_verify['terminal'] or not old_verify['bitwise_metrics'] or sha(audit_path)!=old_verify['original_summary_sha256']:
        raise ValueError('Replayed original physical audit required')
    with np.load(probe/'context.npz') as raw:context=dict(raw)
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    d=RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
                             degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    identity=provenance(['examples.trajectory_spatial_parts_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],
                input_sha256=inputs,parts_directory=str(parts_directory),spec=spec,quadrature=[])
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save();integrals={}
    for order in (16,20):
        sample=TransportSampler(d,mesh,order);integral=lambda value:float(np.sum(sample.weight*value))
        volume=float(sample.weight.sum());zero=np.zeros_like(sample.weight)
        q=dict(order=order,cases=[],weight_sha256=fingerprint(sample.weight));result['quadrature'].append(q);save()
        previous_direction=None;cached={}
        original_q=next(item for item in old_audit['quadrature'] if item['order']==order)
        for case in record['cases']:
            index=case['index'];name=case['direction'];h=case['amplitude']
            if previous_direction!=name:
                cached={};previous_direction=name
                du=context[name+'_fluid'];db=context[name+'_magnetic']
                df,dp,_=expand(sample,context['state'],du,db)
                direction=components(df,dp,sample.electric(context[name+'_current']))
                strong,charge,diagnostics=interior(sample,context['state'],du,context[name+'_current'],
                                                  probe_record['physical_bias'],ref.ion_mass,ref.electron_mass)
                reference=analytic(mode,sample.R,sample.Z,probe_record['directions'][name]['normalization']) if name=='independent_mode' else None
            states={};rates={};original_rates={}
            for sign_row in case['signs']:
                sign=sign_row['sign']
                with np.load(parts_directory/sign_row['arrays_file']) as raw:rates[sign]=dict(raw)
                states[sign]=rates[sign]['state']
                original_sign=next(s for s in probe_record['cases'][index]['signs'] if s['sign']==sign)
                with np.load(probe/original_sign['arrays_file']) as raw:
                    original_rates[sign]={k:raw[k+'_rate' if k!='electric' else k] for k in ('fluid','magnetic','current','electric')}
            full,full_fingerprints=values(sample,states,original_rates,h)
            original_case=next(c for c in original_q['cases'] if c['index']==index)
            if any(value!=original_case['fingerprints'][key] for key,value in full_fingerprints.items()):
                raise ValueError('Original complete physical field fingerprints must match')
            row=dict(index=index,direction=name,amplitude=h,full_original_fingerprints=full_fingerprints,
                     interior_diagnostics=diagnostics,groups=[],sum_vs_full={},amplitude_comparisons=[])
            if reference is not None:
                normalization=probe_record['directions'][name]['normalization'];target=mode.fields(sample.R,sample.Z)
                f0=ref.fields(sample.Z);p0=np.stack((f0['electron_pressure'],f0['ion_pressure']),axis=-1)
                pressure_target=p0*(d.thermal.gammas*target['physical'][...,:1]/f0['number_density'][...,None]+target['physical'][...,4:6])/normalization
                expected=components(mode.growth*target['physical']/normalization,mode.growth*pressure_target,
                                    mode.growth*target['current']/normalization,target['electric']/normalization)
                row['analytic_reference_sum_vs_mode']={key:compare(sample,sum(v.get(key,zero) for v in reference.values()),value)
                                                        for key,value in expected.items()}
            total={key:np.zeros_like(value) for key,value in full.items()}
            for group in GROUPS:
                observed,fingerprints=values(sample,states,rates,h,group+'_')
                metrics={};table=[]
                for key,value in observed.items():
                    squared=integral(value*value);target=direction.get(key)
                    den=integral(target*target) if target is not None else 0.
                    product=integral(value*target) if target is not None else 0.
                    full_den=integral(full[key]*full[key]);full_product=integral(value*full[key])
                    metrics[key]=dict(squared_integral=squared,RMS=float(np.sqrt(squared/volume)),
                        direction_product=product,direction_squared_integral=den,
                        projection_per_second=product/den if den else None,
                        norm_relative_to_full=float(np.sqrt(squared/full_den)) if full_den else None,
                        signed_projection_on_full=full_product/full_den if full_den else None)
                    table.append([squared,den,product,full_den,full_product]);total[key]+=value
                part=dict(name=group,metrics=metrics,fingerprints=fingerprints)
                target=strong.get(group)
                if target is not None:
                    part['strong_interior_errors']={key:compare(sample,observed[key],value) for key,value in target.items()}
                if group in charge:
                    part['charge_volume_alternatives']={kind:{key:compare(sample,observed[key],value)
                        for key,value in target.items() if key in ('electron_entropy','ion_entropy','electron_pressure','ion_pressure')}
                        for kind,target in charge[group].items()}
                if reference is not None:
                    target=reference.get(group,{})
                    part['independent_continuum_errors']={key:compare(sample,value,target.get(key,zero)) for key,value in observed.items()}
                if group in cached:
                    previous_h,previous=cached[group]
                    row['amplitude_comparisons'].append(dict(group=group,first_amplitude=previous_h,second_amplitude=h,
                        actual_field_differences={key:compare(sample,value,previous[key]) for key,value in observed.items()}))
                cached[group]=(h,observed)
                integrals[f'q{order}_case{index:02}_{group}']=np.asarray(table)
                row['groups'].append(part)
            row['sum_vs_full']={key:compare(sample,value,full[key]) for key,value in total.items()}
            row['sum_fingerprints']={key:fingerprint(value) for key,value in total.items()}
            q['cases'].append(row);save();print('physical spatial parts',order,index,name,h,flush=True)
        del cached,reference,strong,charge,direction,total,full,rates,original_rates
    path=destination/'integrals.npz';np.savez_compressed(path,**integrals)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.trajectory_spatial_parts_audit'])!=identity:
        raise ValueError('Original physical audit inputs/source changed')
    result.update(terminal=True,all_cases_complete=True,arrays_file=path.name,arrays_sha256=sha(path),
                  scalar_names=[*NAMES,'electron_pressure','ion_pressure',
                    *['current_'+axis for axis in ('radial','azimuthal','axial')],
                    *['electric_'+axis for axis in ('radial','azimuthal','axial')]],
                  integral_columns=['squared','direction_squared','direction_product','full_squared','full_product'],
                  scope='Separate physical contributions and independent continuum control; strong interior comparisons omit face distributions. E is algebraic. No causal, eigenvalue or formation claim.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('parts_directory')
    a=p.parse_args();run(a.destination,a.parts_directory)
