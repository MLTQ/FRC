"""Full-vector entropy interface and potential-transfer diagnostic."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_species import RegularSpecies
from .continuum_contact_mode import ContinuumContactMode
from .species_weak_sampling import WeakSampler
from .trajectory_current_charge import face_edges
from .trajectory_tangent_probe import verified_inputs
from .trajectory_tangent_audit import compare
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def relative(value,target):
    den=float(np.linalg.norm(target));difference=value-target
    return dict(relative=float(np.linalg.norm(difference)/den) if den else None,
                max_absolute=float(np.max(abs(difference))))


def analytic_control(weak):
    s=weak.sample;d=weak.fluid;k=2*np.pi/d.length
    entropy=(d.r/d.radius)**2+.2*np.cos(k*d.z)
    exact=(s.R/d.radius)**2+.2*np.cos(k*s.Z)
    radial=np.broadcast_to(s.R*k*np.cos(k*s.Z),exact.shape)
    axial=np.broadcast_to(-2*np.sin(k*s.Z),exact.shape)
    j=np.stack((radial,np.zeros_like(radial),axial),axis=-1)
    jr=weak.rf[:,None]*k*np.cos(k*s.z[None,:])
    jz=np.broadcast_to(-2*np.sin(k*weak.zf)[None,:],(len(s.r),len(weak.zf)))
    terms=weak.decompose(entropy,j,np.zeros_like(exact),jr,jz,weak.projection)
    strong=2*s.R**2*k*np.cos(k*s.Z)/d.radius**2+.4*k*np.sin(k*s.Z)**2
    target=weak.adjoint(strong)
    return terms,dict(identity=relative(terms['sum'],terms['direct']),
        scalar=relative(weak.scalar(entropy),exact),
        radial_gradient=relative(weak.DR@weak.shape(entropy)@weak.Z.T,np.broadcast_to(2*s.R/d.radius**2,exact.shape)),
        axial_gradient=relative(weak.R@weak.shape(entropy)@weak.DZ.T,np.broadcast_to(-.2*k*np.sin(k*s.Z),exact.shape)),
        advective_vector=relative(terms['advection'],target),
        analytic_target_fingerprint=fingerprint(strong))


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh species-interface output required')
    parts=Path('outputs/trajectory_spatial_parts');probe=Path('outputs/trajectory_tangent_probe')
    record,inputs=verified_inputs(parts,str(parts)+'_replay_verification')
    charge=Path('outputs/trajectory_current_charge')
    charge_record,charge_inputs=verified_inputs(charge,str(charge)+'_replay_verification')
    inputs.update(charge_inputs);inputs.update(record['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Unchanged replayed inputs required')
    probe_record=json.loads((probe/'summary.json').read_text())
    with np.load(probe/'context.npz') as raw:context=dict(raw)
    if sha(probe/'context.npz')!=probe_record['context_sha256']:raise ValueError('Original direction context required')
    inputs[str(probe/'context.npz')]=probe_record['context_sha256']
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    d=RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
                             degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    base=context['state'];geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    model=RegularCurrent(RegularFieldMap(d,geometry),base[...,0],ref.ion_mass)
    species=RegularSpecies(model);production=(species.potential.radial,species.potential.axial)
    identity=provenance(['examples.species_interface_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,
                spec=spec,production_order=12,quadrature=[],identity_reporting_threshold=1e-9)
    def save():
        path=destination/'summary.pending.json';path.write_text(json.dumps(result,indent=2)+'\n');path.replace(destination/'summary.json')
    save();identity_errors=[]
    count=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    ex=face_edges(.5*mesh.rf**2);ez=face_edges(mesh.zf,True)
    for order in (12,16,20):
        weak=WeakSampler(d,mesh,order);sample=weak.sample
        q=dict(order=order,cases=[],projection_vs_production={axis:relative(value,target)
            for axis,value,target in zip(('radial','axial'),weak.projection,production)})
        result['quadrature'].append(q)
        manufactured,q['analytic_control']=analytic_control(weak)
        path=destination/f'control_q{order}.npz'
        np.savez_compressed(path,**manufactured,projection_radial=weak.projection[0],projection_axial=weak.projection[1])
        q.update(control_arrays_file=path.name,control_arrays_sha256=sha(path));save()
        n=sample.scalar(base[...,0])
        for name in ('saved_error','independent_mode'):
            current=context[name+'_current'];j,gradient=sample.electric(current,gradient=True)
            div=np.trace(gradient,axis1=-2,axis2=-1)
            jr=weak.rf[:,None]*(ex@current[:count].reshape(mesh.nr,mesh.nz)@sample.nz.T)
            jz=sample.nx@current[count+theta:].reshape(mesh.nr+1,mesh.nz)@ez.T
            row=dict(direction=name,current_fingerprint=fingerprint(j),divergence_fingerprint=fingerprint(div),species=[])
            if order in (16,20):
                previous=next(c for c in charge_record['cases'] if c['order']==order and c['direction']==name)
                with np.load(charge/previous['arrays_file']) as raw:
                    exact=(np.array_equal(div,raw['divergence']) and
                           np.array_equal(jr,raw['radial_traces'].reshape(jr.shape)) and
                           np.array_equal(jz,raw['axial_traces'].reshape(jz.shape)))
                if not exact or fingerprint(j)!=previous['current_fingerprint']:raise ValueError('Original charge reconstruction required')
                row['previous_charge_arrays_exact']=True
            raw_vectors={};direction=context[name+'_fluid']
            for index,(label,a) in enumerate(zip(('electron','ion'),(-model.beta,model.alpha))):
                slot=index+4;entropy=base[...,slot]/base[...,0]
                original=species.transport(current,entropy).reshape(-1,d.nz)
                multiplier=-a/E/d.volume.reshape(-1,d.nz)
                original_rate=(multiplier*original).reshape(d.shape[:-1])
                physical=sample.scalar(original_rate)/n
                quotient=sample.scalar(base[...,slot])/n
                specific_direction=(sample.scalar(direction[...,slot])-quotient*sample.scalar(direction[...,0]))/n
                norm=float(np.sum(sample.weight*specific_direction**2))
                item=dict(name=label,representation={},projections={},amplitude_comparisons=[])
                pair=weak.shape(entropy)*original
                item['original_native_conservation']=dict(
                    sum_K=float(original.sum()),sum_absolute_K=float(np.sum(abs(original))),
                    entropy_K_product=float(pair.sum()),sum_absolute_entropy_K=float(np.sum(abs(pair))),
                    entropy_density_rate_integral=float(np.sum(d.volume*original_rate)),
                    physical_entropy_density_rate_integral=float(np.sum(sample.weight*sample.scalar(original_rate))))
                raw_vectors[label+'_production_K']=original
                raw_vectors[label+'_production_Cdot']=original_rate
                for axis,derivative in [('value',None),('radial',0),('axial',2)]:
                    if derivative is None:
                        nodal=sample.scalar(entropy);other=quotient
                    else:
                        nodal=sample.scalar(entropy,axis=derivative)
                        other=(sample.scalar(base[...,slot],axis=derivative)-quotient*sample.scalar(base[...,0],axis=derivative))/n
                    item['representation'][axis]=compare(sample,nodal,other)
                nodal_strong=-a/E*(j[...,0]*sample.scalar(entropy,axis=0)+j[...,2]*sample.scalar(entropy,axis=2))/n
                quotient_strong=-a/E*sum(j[...,axis]*(sample.scalar(base[...,slot],axis=axis)-quotient*sample.scalar(base[...,0],axis=axis))/n for axis in (0,2))/n
                item['strong_comparisons']={'nodal_entropy':compare(sample,physical,nodal_strong),
                                            'quotient_entropy':compare(sample,physical,quotient_strong)}
                for projection_label,projection in [('production',production),('independent',weak.projection)]:
                    terms=weak.decompose(entropy,j,div,jr,jz,projection)
                    identity_check=relative(terms['sum'],terms['direct'])
                    identity_errors.append(identity_check['relative'])
                    report=dict(identity=identity_check,direct_vs_original=relative(terms['direct'],original),components={})
                    for term,value in terms.items():
                        raw_vectors[f'{label}_{projection_label}_{term}']=value
                        rate=sample.scalar((multiplier*value).reshape(d.shape[:-1]))/n
                        report['components'][term]=dict(squared_integral=float(np.sum(sample.weight*rate**2)),
                            direction_product=float(np.sum(sample.weight*rate*specific_direction)),
                            direction_squared=norm,projection_per_second=float(np.sum(sample.weight*rate*specific_direction)/norm) if norm else None,
                            original_comparison=compare(sample,rate,physical),fingerprint=fingerprint(rate))
                    item['projections'][projection_label]=report
                for case in record['cases']:
                    if case['direction']!=name:continue
                    states={};rates={}
                    for sign in ('minus','plus'):
                        with np.load(parts/f"case_{case['index']:02}_{sign}.npz") as data:
                            states[sign]=data['state'];rates[sign]=data[label+'_species_fluid'][...,slot]
                    h=case['amplitude'];centered=(rates['plus']-rates['minus'])/(2*h)
                    actual=(sample.scalar(rates['plus'])/sample.scalar(states['plus'][...,0])-
                            sample.scalar(rates['minus'])/sample.scalar(states['minus'][...,0]))/(2*h)
                    item['amplitude_comparisons'].append(dict(amplitude=h,
                        native=relative(centered,original_rate),physical=compare(sample,actual,physical)))
                    raw_vectors[f'{label}_centered_h{h}']=centered
                row['species'].append(item)
            path=destination/f'{name}_q{order}.npz';np.savez_compressed(path,**raw_vectors)
            row.update(arrays_file=path.name,arrays_sha256=sha(path));q['cases'].append(row);save()
            print('whole-vector weak entropy',order,name,'maximum identity error',max(identity_errors),flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.species_interface_audit'])!=identity:
        raise ValueError('Audit source or input changed')
    result.update(terminal=True,all_cases_complete=True,maximum_identity_relative=max(identity_errors),
                  all_identity_comparisons_below_reporting_threshold=max(identity_errors)<1e-9,
                  scope='Complete stationary current-entropy weak budget, not a replacement interface flux or a causal growth proof.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
