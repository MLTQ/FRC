"""Independent physical fields, continuum control and tangent-amplitude audit."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .continuum_wall_probe import physical
from .current_transport_sampling import TransportSampler
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha
from .trajectory_tangent_probe import verified_inputs
from .coupled_rate_physical import response,power
from .split_mode_audit import fixture
from plasma_bench.regular_current import RegularCurrent


NAMES=['density','radial_velocity','azimuthal_velocity','axial_velocity',
    'electron_entropy','ion_entropy','radial_B','azimuthal_B','axial_B']


def expand(sample,state,fluid,magnetic):
    fields=physical(sample,state,SimpleNamespace(fluid=fluid,magnetic=magnetic))
    n=sample.scalar(state[...,0])
    entropy=np.stack([sample.scalar(state[...,k])/n for k in (4,5)],axis=-1)
    pressure=E*n[...,None]*np.exp(entropy)*(n[...,None]/sample.d.thermal.n0)**(sample.d.thermal.gammas-1)
    pressure_rate=pressure*(sample.d.thermal.gammas*fields[...,:1]/n[...,None]+fields[...,4:6])
    if not np.isfinite(fields).all() or not np.isfinite(pressure_rate).all() or np.any(pressure<=0):
        raise ValueError('Finite physical rates and positive sampled pressures required')
    return fields,pressure_rate,dict(density_min=float(n.min()),
        temperature_min_by_species=np.min(pressure/(E*n[...,None]),axis=(0,1)).tolist())


def components(fields,pressure,current,electric=None):
    result={name:fields[...,k] for k,name in enumerate(NAMES)}
    result.update(electron_pressure=pressure[...,0],ion_pressure=pressure[...,1])
    result.update({f'current_{axis}':current[...,k] for k,axis in enumerate(('radial','azimuthal','axial'))})
    if electric is not None:result.update({f'electric_{axis}':electric[...,k] for k,axis in enumerate(('radial','azimuthal','axial'))})
    return result


def compare(sample,value,target):
    difference=value-target
    def ratio(mask):
        den=float(np.sum(sample.weight[mask]*target[mask]**2))
        num=float(np.sum(sample.weight[mask]*difference[mask]**2))
        return float(np.sqrt(num/den)) if den else None
    return dict(relative_L2=ratio(slice(None)),max_absolute=float(np.max(abs(difference))),
        axis_element=ratio(sample.r<=sample.d.radius/sample.d.elements),
        wall_element=ratio(sample.r>=sample.d.radius*(1-1/sample.d.elements)))


def run(destination,probe_directory):
    destination=Path(destination);probe_directory=Path(probe_directory)
    if destination.exists():raise ValueError('Fresh tangent-audit destination required')
    probe=json.loads((probe_directory/'summary.json').read_text())
    source=json.loads((probe_directory/'provenance.json').read_text())
    if not probe['terminal'] or 'exception_type' in probe:raise ValueError('Complete original tangent protocol required')
    if any(sha(p)!=h for p,h in source['source_files'].items()):raise ValueError('Original tangent source required')
    _,inputs=verified_inputs(probe_directory,str(probe_directory)+'_replay_verification')
    inputs.update(probe['input_sha256'])
    def load(name,digest):
        p=probe_directory/name
        if sha(p)!=digest:raise ValueError('Original complete tangent raw bytes required')
        inputs[str(p)]=digest
        with np.load(p) as raw:return {k:raw[k].copy() for k in raw.files}
    context=load(probe['context_file'],probe['context_sha256'])
    base_rate=load(probe['base']['arrays_file'],probe['base']['arrays_sha256'])
    loaded=[];rejected=[]
    for index,case in enumerate(probe['cases']):
        states=load(case['input_file'],case['input_sha256'])
        if case['status']!='complete':rejected.append(dict(index=index,record=case));continue
        rates={r['sign']:load(r['arrays_file'],r['arrays_sha256']) for r in case['signs']}
        action=load(case['action_file'],case['action_sha256']);loaded.append((index,case,states,rates,action))
    baseline_directory=Path('outputs/trajectory_tangent_probe')
    baseline,baseline_inputs=verified_inputs(baseline_directory,str(baseline_directory)+'_replay_verification')
    inputs.update(baseline_inputs)
    baseline_loaded={}
    for index,case in enumerate(baseline['cases']):
        with np.load(baseline_directory/case['input_file']) as raw:states=dict(raw)
        rates={}
        for sign in case['signs']:
            with np.load(baseline_directory/sign['arrays_file']) as raw:rates[sign['sign']]=dict(raw)
        baseline_loaded[index]=(states,rates)
    old_audit_path=Path('outputs/trajectory_tangent_audit/summary.json')
    old_audit=json.loads(old_audit_path.read_text())
    old_verifier=Path('outputs/trajectory_tangent_audit_replay_verification/summary.json')
    verifier=json.loads(old_verifier.read_text())
    if not verifier['terminal'] or not verifier['bitwise_metrics'] or sha(old_audit_path)!=verifier['original_summary_sha256']:
        raise ValueError('Complete original physical tangent audit required')
    inputs[str(old_audit_path)]=sha(old_audit_path);inputs[str(old_verifier)]=sha(old_verifier)
    mode=ContinuumContactMode();ref=mode.ref;spec=probe['spec']
    if mode.input_sha256!=probe['reference_input_sha256']:raise ValueError('Original independent continuum mode required')
    d=RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
        degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    coarse=RegularCoordinateFluid(ref.radius,ref.length,8,spec['points'],ref.ion_mass+ref.electron_mass,
        degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    _,_,production_base,system,bias,_,_=fixture(spec)
    if not np.array_equal(production_base,context['state']):raise ValueError('Exact original base required for current reconstruction')
    model=RegularCurrent(system.mapping,production_base[...,0],ref.ion_mass)
    input_currents={}
    for index,case,states,rates,action in loaded:
        for sign in ('minus','plus'):input_currents[f'case_{index}_{sign}']=model.current(states[sign+'_magnetic'])
    input_currents['base']=model.current(context['magnetic'])
    identity=provenance(['examples.coupled_fluid_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    current_path=destination/'input_currents.npz';np.savez_compressed(current_path,**input_currents)
    result=dict(terminal=False,source_sha256=identity['source_sha256'],probe_directory=str(probe_directory),
        input_sha256=inputs,spec=spec,all_cases_complete=probe['all_cases_complete'],rejected_cases=rejected,quadrature=[],
        input_currents_file=current_path.name,input_currents_sha256=sha(current_path))
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save();raw_integrals={}
    for order in (16,20):
        sample=TransportSampler(d,mesh,order);coarse_sample=TransportSampler(coarse,mesh,order);target=mode.fields(sample.R,sample.Z)
        original_q=next(q for q in old_audit['quadrature'] if q['order']==order)
        background=ref.fields(sample.Z);p0=np.stack((background['electron_pressure'],background['ion_pressure']),axis=-1)
        target_pressure=p0*(d.thermal.gammas*target['physical'][...,:1]/background['number_density'][...,None]+target['physical'][...,4:6])
        base_fields,base_pressure,base_positivity=expand(sample,context['state'],base_rate['fluid_rate'],base_rate['magnetic_rate'])
        base_values=components(base_fields,base_pressure,sample.electric(base_rate['current_rate']),sample.electric(base_rate['electric']))
        integral=lambda value:float(np.sum(sample.weight*value))
        q=dict(order=order,cases=[],amplitude_comparisons=[],base_positivity=base_positivity,
            base_physical_RMS={k:float(np.sqrt(integral(v*v)/sample.weight.sum())) for k,v in base_values.items()},
            base_fingerprints={k:fingerprint(v) for k,v in base_values.items()},weight_sha256=fingerprint(sample.weight))
        q['base_physical_conservation']=power(sample,context['state'],context['magnetic'],input_currents['base'],base_rate,bias.axial_bias,ref.ion_mass,ref.electron_mass)
        cached={}
        for index,case,states,rates,action in loaded:
            name=case['direction'];amplitude=case['amplitude'];normalization=probe['directions'][name]['normalization']
            du=context[name+'_fluid'];db=context[name+'_magnetic']
            direction_fields,direction_pressure,_=expand(sample,context['state'],du,db)
            direction_current=sample.electric(context[name+'_current'])
            direction=components(direction_fields,direction_pressure,direction_current)
            pushed_fields,pushed_pressure,_=expand(sample,context['state'],action['fluid'],action['magnetic'])
            pushed_current=sample.electric(action['current']);pushed=components(pushed_fields,pushed_pressure,pushed_current)
            minus=expand(sample,states['minus_state'],rates['minus']['fluid_rate'],rates['minus']['magnetic_rate'])
            plus=expand(sample,states['plus_state'],rates['plus']['fluid_rate'],rates['plus']['magnetic_rate'])
            rhs_fields=(plus[0]-minus[0])/(2*amplitude);rhs_pressure=(plus[1]-minus[1])/(2*amplitude)
            current=sample.electric((rates['plus']['current_rate']-rates['minus']['current_rate'])/(2*amplitude))
            electric=sample.electric((rates['plus']['electric']-rates['minus']['electric'])/(2*amplitude))
            values=components(rhs_fields,rhs_pressure,current,electric);cached.setdefault(name,[]).append((amplitude,values))
            alignment={};squared={};products={}
            for key,value in values.items():
                energy=integral(value*value);squared[key]=energy
                if key in direction:
                    tangent=direction[key];den=integral(tangent*tangent);product=integral(value*tangent);products[key]=product
                    projection=product/den if den else None
                    residual=integral((value-projection*tangent)**2) if den else None
                    alignment[key]=dict(direction_squared_integral=den,rhs_squared_integral=energy,
                        projection_per_second=projection,orthogonal_relative=float(np.sqrt(residual/energy)) if den and energy else None,
                        residual_to_coordinate_projection=compare(sample,value,case['coordinate_projection']*tangent))
            row=dict(index=index,direction=name,amplitude=amplitude,alignment=alignment,
                physical_rhs_squared_integrals=squared,positivity=[minus[2],plus[2]],
                physical_rhs_vs_conservative_pushforward={k:compare(sample,values[k],v) for k,v in pushed.items()},
                fingerprints={k:fingerprint(v) for k,v in dict(direction_fields=direction_fields,direction_pressure=direction_pressure,
                    direction_current=direction_current,pushed_fields=pushed_fields,pushed_pressure=pushed_pressure,pushed_current=pushed_current,
                    rhs_fields=rhs_fields,rhs_pressure=rhs_pressure,current_rate=current,electric=electric).items()})
            coarse_states,coarse_rates=baseline_loaded[index]
            baseline_values,baseline_fingerprints=response(coarse_sample,coarse_states,coarse_rates,amplitude)
            original_case=next(c for c in original_q['cases'] if c['index']==index)
            if any(original_case['fingerprints'][k]!=v for k,v in baseline_fingerprints.items()):
                raise ValueError('Original full physical baseline fingerprints must match exactly')
            row['original_F8_fingerprints_exact']=True
            row['physical_rhs_vs_F8']={k:compare(sample,value,baseline_values[k]) for k,value in values.items()}
            row['mesh_difference_fingerprints']={k:fingerprint(value-baseline_values[k]) for k,value in values.items()}
            row['instantaneous_conservation']={}
            for sign in ('minus','plus'):
                physical=power(sample,states[sign+'_state'],states[sign+'_magnetic'],input_currents[f'case_{index}_{sign}'],
                               rates[sign],bias.axial_bias,ref.ion_mass,ref.electron_mass)
                native=next(r for r in case['signs'] if r['sign']==sign)['diagnostics']
                physical['native_powers_W']=native['powers'];physical['native_total_power_W']=native['total_power']
                row['instantaneous_conservation'][sign]=physical
            if name=='independent_mode':
                independent=components(mode.growth*target['physical']/normalization,
                    mode.growth*target_pressure/normalization,mode.growth*target['current']/normalization,target['electric']/normalization)
                initial_target=components(target['physical']/normalization,target_pressure/normalization,target['current']/normalization)
                row.update(continuum_rate_errors={k:compare(sample,values[k],v) for k,v in independent.items()},
                    initial_representation_errors={k:compare(sample,direction[k],v) for k,v in initial_target.items()},
                    continuum_fingerprints={k:fingerprint(v) for k,v in independent.items()})
            q['cases'].append(row)
            raw_integrals[f'q{order}_case{index}_rhs_squared']=np.array(list(squared.values()))
            raw_integrals[f'q{order}_case{index}_direction_squared']=np.array([v['direction_squared_integral'] for v in alignment.values()])
            raw_integrals[f'q{order}_case{index}_rhs_direction_products']=np.array(list(products.values()))
        for name,observed in cached.items():
            for (first,a),(second,b) in zip(observed[:-1],observed[1:]):
                q['amplitude_comparisons'].append(dict(direction=name,first_amplitude=first,second_amplitude=second,
                    differences_relative_to_first={k:compare(sample,b[k],v) for k,v in a.items()}))
        result['quadrature'].append(q);save();print('independent physical tangent audit',order,flush=True)
    path=destination/'integrals.npz';np.savez_compressed(path,**raw_integrals)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.coupled_fluid_audit'])!=identity:
        raise ValueError('Physical audit input or source changed')
    result.update(terminal=True,arrays_file=path.name,arrays_sha256=sha(path),
        scalar_rhs_names=[*NAMES,'electron_pressure','ion_pressure',
            *['current_'+axis for axis in ('radial','azimuthal','axial')],*['electric_'+axis for axis in ('radial','azimuthal','axial')]],
        scope='Independent complete physical reconstruction, F8 response comparison and instantaneous original/physical conservation for fluid-space refinement. Physical-RHS derivatives differ from conservative-Jacobian pushforwards through finite amplitude and base stationarity defect. Electric response is algebraic, not an evolved electric time derivative. Component projections are directional diagnostics, not eigenvalue or causal claims. No unit tests.')
    save();return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('probe_directory')
    args=parser.parse_args();run(args.destination,args.probe_directory)
