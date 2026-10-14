"""Actual full-response field and conservation comparisons on three fluid grids."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.regular_current import RegularCurrent
from .current_transport_sampling import TransportSampler
from .split_mode_audit import fixture
from .trajectory_tangent_probe import verified_inputs
from .trajectory_tangent_audit import compare
from .coupled_rate_physical import response,power
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh full-grid comparison output required')
    datasets={};inputs={}
    for elements in (8,16,32):
        probe_dir=Path('outputs')/('trajectory_tangent_probe' if elements==8 else f'coupled_fluid_probe_F{elements}')
        audit_dir=Path('outputs')/('trajectory_tangent_audit' if elements==8 else f'coupled_fluid_audit_F{elements}')
        probe,checked=verified_inputs(probe_dir,str(probe_dir)+'_replay_verification');inputs.update(checked);inputs.update(probe['input_sha256'])
        audit,checked=verified_inputs(audit_dir,str(audit_dir)+'_replay_verification');inputs.update(checked);inputs.update(audit['input_sha256'])
        if not probe['all_cases_complete'] or not audit['all_cases_complete']:raise ValueError('All declared coupled cases must be present')
        with np.load(probe_dir/probe['context_file']) as raw:context=dict(raw)
        archived_currents={}
        if elements!=8:
            with np.load(audit_dir/audit['input_currents_file']) as raw:archived_currents=dict(raw)
        cases=[]
        for item in probe['cases']:
            with np.load(probe_dir/item['input_file']) as raw:states=dict(raw)
            rates={}
            for sign in item['signs']:
                with np.load(probe_dir/sign['arrays_file']) as raw:rates[sign['sign']]=dict(raw)
            cases.append((item,states,rates))
        datasets[elements]=dict(probe=probe,audit=audit,context=context,cases=cases,archived_currents=archived_currents)
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original complete comparison inputs required')
    mode,d,base,system,bias,_,_=fixture(datasets[8]['probe']['spec']);ref=mode.ref
    model=RegularCurrent(system.mapping,base[...,0],ref.ion_mass);mesh=system.geometry.mesh
    identity=provenance(['examples.coupled_refinement_comparison']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,quadrature=[])
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save();integrals={}
    for order in (16,20):
        samples={elements:TransportSampler(RegularCoordinateFluid(ref.radius,ref.length,elements,64,
            ref.ion_mass+ref.electron_mass,degree=5,gammas=(ref.gamma_e,ref.gamma_i)),mesh,order) for elements in datasets}
        sample=samples[8];q=dict(order=order,cases=[]);result['quadrature'].append(q)
        for index in range(6):
            values={};row=dict(index=index,grids={},pairs=[])
            for elements,data in datasets.items():
                case,states,rates=data['cases'][index]
                if elements==8:row.update(direction=case['direction'],amplitude=case['amplitude'])
                elif (case['direction'],case['amplitude'])!=(row['direction'],row['amplitude']):raise ValueError('Same fixed direction and amplitude required')
                observed,prints=response(samples[elements],states,rates,case['amplitude']);values[elements]=observed
                oldq=next(x for x in data['audit']['quadrature'] if x['order']==order)
                previous=next(c for c in oldq['cases'] if c['index']==index)
                if any(previous['fingerprints'][k]!=v for k,v in prints.items()):raise ValueError('All complete physical response fingerprints must match')
                current_record=dict(response_fingerprints_exact=True,conservation={})
                for sign in ('minus','plus'):
                    current=model.current(states[sign+'_magnetic'])
                    if elements!=8 and not np.array_equal(current,data['archived_currents'][f'case_{index}_{sign}']):raise ValueError('Exact signed input current required')
                    conservation=power(samples[elements],states[sign+'_state'],states[sign+'_magnetic'],current,rates[sign],bias.axial_bias,ref.ion_mass,ref.electron_mass)
                    if elements!=8 and any(previous['instantaneous_conservation'][sign][k]!=v for k,v in conservation.items()):
                        raise ValueError('Every original physical conservation metric must replay exactly')
                    native=next(s for s in case['signs'] if s['sign']==sign)['diagnostics']
                    conservation.update(native_powers_W=native['powers'],native_total_power_W=native['total_power'])
                    current_record['conservation'][sign]=conservation
                row['grids'][str(elements)]=current_record
            for first,second in [(8,16),(16,32),(8,32)]:
                differences={};table=[];prints={}
                for key,target in values[first].items():
                    value=values[second][key];delta=value-target
                    differences[key]=compare(sample,value,target);prints[key]=fingerprint(delta)
                    table.append([np.sum(sample.weight*delta**2),np.sum(sample.weight*target**2),
                                  np.sum(sample.weight*value**2),np.sum(sample.weight*value*target)])
                key=f'q{order}_case{index}_F{first}_F{second}';integrals[key]=np.asarray(table)
                row['pairs'].append(dict(first_elements=first,second_elements=second,differences=differences,
                                         difference_fingerprints=prints,integrals_key=key))
            q['cases'].append(row);save();print('full physical mesh comparison',order,index,row['direction'],row['amplitude'],flush=True)
    path=destination/'integrals.npz';np.savez_compressed(path,**integrals)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.coupled_refinement_comparison'])!=identity:raise ValueError('Comparison source/input changed')
    result.update(terminal=True,all_cases_complete=True,arrays_file=path.name,arrays_sha256=sha(path),
        scalar_names=list(values[8]),integral_columns=['difference_squared','coarse_squared','fine_squared','signed_product'],
        scope='Full fixed-direction mesh response and instantaneous conservation; not coupled eigenvalue, trajectory, local-flux or formation acceptance.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
