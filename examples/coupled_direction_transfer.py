"""Physical audit before evaluating full RHS on nested fluid grids."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler
from .subdivide_fluid_direction import transfer,parent_faces
from .trajectory_tangent_probe import verified_inputs
from .trajectory_tangent_audit import expand,components,compare
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def finite(sample,state,magnetic,current,bias):
    n,v,_=sample.fluid(state)
    entropy=np.stack([sample.scalar(state[...,k])/n for k in (4,5)],axis=-1)
    p=E*n[...,None]*np.exp(entropy)*(n[...,None]/sample.d.thermal.n0)**(sample.d.thermal.gammas-1)
    b=sample.magnetic(magnetic);b[...,2]+=bias
    fields=np.concatenate((n[...,None],v,entropy,b),axis=-1)
    return components(fields,p,sample.electric(current))


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh direction-transfer output required')
    probe=Path('outputs/trajectory_tangent_probe');record,inputs=verified_inputs(probe,str(probe)+'_replay_verification')
    inputs.update(record['input_sha256'])
    context_path=probe/record['context_file']
    if sha(context_path)!=record['context_sha256']:raise ValueError('Original complete direction context required')
    inputs[str(context_path)]=sha(context_path)
    with np.load(context_path) as raw:context=dict(raw)
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    def fluid(elements):
        return RegularCoordinateFluid(ref.radius,ref.length,elements,spec['points'],ref.ion_mass+ref.electron_mass,
                                     degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    coarse=fluid(8);mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    base=context['state']
    if not np.array_equal(base,np.broadcast_to(base[0,0],base.shape)):raise ValueError('Exactly radially constant base required')
    identity=provenance(['examples.coupled_direction_transfer']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,
                original_spec=spec,directions=record['directions'],amplitudes=record['amplitudes'],physical_bias=record['physical_bias'],grids=[],transfer_relative_bound=1e-9)
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save();maximum=0.
    for elements in (16,32):
        fine=fluid(elements);fine_base=np.broadcast_to(base[0,0],fine.shape).copy()
        arrays=dict(state=fine_base,magnetic=context['magnetic']);grid=dict(elements=elements,quadrature=[],parent_traces={},finite_transfer=[])
        for name in ('saved_error','independent_mode'):
            du=transfer(coarse,fine,context[name+'_fluid'])
            for key,value in [('fluid',du),('magnetic',context[name+'_magnetic']),('current',context[name+'_current'])]:arrays[name+'_'+key]=value
            before=parent_faces(context[name+'_fluid'],8);after=parent_faces(du,8)
            if not np.array_equal(before,after):raise ValueError('Exact owning-parent conservative face data required')
            if np.any(du[0,0,:,1:3]) or np.any(du[-1,-1,:,1]):raise ValueError('Exact original inactive/wall constraints required')
            grid['parent_traces'][name]=dict(bitwise_equal=True,fingerprint=fingerprint(after))
        path=destination/f'F{elements}_context.npz';np.savez_compressed(path,**arrays)
        grid.update(context_file=path.name,context_sha256=sha(path));result['grids'].append(grid);save()
        for order in (16,20):
            old=TransportSampler(coarse,mesh,order);new=TransportSampler(fine,mesh,order)
            q=dict(order=order,directions=[])
            for name in ('saved_error','independent_mode'):
                du=arrays[name+'_fluid'];db=arrays[name+'_magnetic'];current=arrays[name+'_current']
                before=expand(old,base,context[name+'_fluid'],db);after=expand(new,fine_base,du,db)
                old_values=components(before[0],before[1],old.electric(current));new_values=components(after[0],after[1],new.electric(current))
                differences={k:compare(new,new_values[k],v) for k,v in old_values.items()}
                for k,v in differences.items():
                    if v['relative_L2'] is not None:maximum=max(maximum,v['relative_L2'])
                    elif v['max_absolute']!=0:raise ValueError('Zero differential component did not transfer exactly: '+k)
                if maximum>result['transfer_relative_bound']:raise ValueError('Physical differential transfer exceeds declared bound')
                if not np.array_equal(old.electric(current),new.electric(current)) or not np.array_equal(old.magnetic(db),new.magnetic(db)):
                    raise ValueError('Bitwise unchanged magnetic/current physical fields required')
                row=dict(direction=name,differences=differences,old_fingerprints={k:fingerprint(v) for k,v in old_values.items()},
                         new_fingerprints={k:fingerprint(v) for k,v in new_values.items()},finite_states=[])
                for h in record['amplitudes']:
                    for sign in (-1,1):
                        lhs=finite(old,base+sign*h*context[name+'_fluid'],sign*h*db,sign*h*current,record['physical_bias'])
                        rhs=finite(new,fine_base+sign*h*du,sign*h*db,sign*h*current,record['physical_bias'])
                        metrics={}
                        for key,value in rhs.items():
                            delta=value-lhs[key];reference=h*old_values[key]
                            den=float(np.sum(new.weight*reference**2));num=float(np.sum(new.weight*delta**2))
                            metrics[key]=dict(relative_to_linear_perturbation=float(np.sqrt(num/den)) if den else None,max_absolute=float(np.max(abs(delta))))
                        row['finite_states'].append(dict(amplitude=h,sign=sign,differences=metrics))
                q['directions'].append(row)
            grid['quadrature'].append(q);save();print('complete physical transfer',elements,order,'maximum relative',maximum,flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.coupled_direction_transfer'])!=identity:raise ValueError('Transfer source/input changed')
    result.update(terminal=True,all_cases_complete=True,maximum_differential_relative=maximum,
                  scope='Conservative diagnostic direction subdivision; finite primitive interpolation differences remain measured. No evolved state or RHS change.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
