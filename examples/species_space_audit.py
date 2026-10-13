"""Fixed-current entropy-source refinement across compatible fluid spaces."""
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
from .species_interface_audit import relative
from .species_mass_audit import mass_reference
from .trajectory_tangent_probe import verified_inputs
from .trajectory_tangent_audit import compare
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh entropy-space audit required')
    prior=Path('outputs/species_mass_audit');record,inputs=verified_inputs(prior,str(prior)+'_replay_verification')
    inputs.update(record['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Unchanged replayed mass audit required')
    with np.load('outputs/trajectory_tangent_probe/context.npz') as raw:context=dict(raw)
    original_base=context['state']
    if not np.array_equal(original_base,np.broadcast_to(original_base[0,0],original_base.shape)):
        raise ValueError('Original base must be exactly radially constant for this transfer')
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    def fluid(elements):
        return RegularCoordinateFluid(ref.radius,ref.length,elements,spec['points'],ref.ion_mass+ref.electron_mass,
                                      degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    original_fluid=fluid(8)
    identity=provenance(['examples.species_space_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,
                original_spec=spec,fluid_elements=[8,16,32],field_mesh_unchanged=True,base_replication_exact=True,cases=[])
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save()
    for order in (16,20):
        original_sampler=WeakSampler(original_fluid,mesh,order).sample
        original_n=original_sampler.scalar(original_base[...,0])
        for elements in result['fluid_elements']:
            d=fluid(elements);base=np.broadcast_to(original_base[0,0],d.shape).copy()
            model=RegularCurrent(RegularFieldMap(d,geometry),base[...,0],ref.ion_mass);species=RegularSpecies(model)
            weak=WeakSampler(d,mesh,order);sample=weak.sample;n=sample.scalar(base[...,0])
            mr=weak.R.T@(weak.wr[:,None]*weak.R);mz=weak.Z.T@(weak.wz[:,None]*weak.Z)
            er,ez,_,_,vr=mass_reference(d)
            for name in ('saved_error','independent_mode'):
                current=context[name+'_current'];j=sample.electric(current)
                if not np.array_equal(j,original_sampler.electric(current)):raise ValueError('Exact unchanged physical current required')
                with np.load(prior/f'{name}_q{order}.npz') as archive:previous=dict(archive)
                row=dict(order=order,elements=elements,direction=name,current_fingerprint=fingerprint(j),
                    base_density_vs_original=compare(sample,n,original_n),
                    radial_mass_analytic=relative(mr/np.sqrt(vr[:,None]*vr[None,:]),er),
                    axial_mass_analytic=relative(mz/d.dz,ez),species=[])
                arrays=dict(current_coefficients=current,base=base,radial_mass=mr,axial_mass=mz)
                for slot,label,a in [(4,'electron',-model.beta),(5,'ion',model.alpha)]:
                    entropy=base[...,slot]/base[...,0];original_entropy=original_base[...,slot]/original_base[...,0]
                    fixed_source=-a/E*sum(j[...,axis]*original_sampler.scalar(original_entropy,axis=axis) for axis in (0,2))
                    reconstructed_source=-a/E*sum(j[...,axis]*sample.scalar(entropy,axis=axis) for axis in (0,2))
                    target=fixed_source/original_n
                    original=((-a/E/d.volume)*species.transport(current,entropy)).reshape(-1,d.nz)
                    adv=-a/E*weak.adjoint(j[...,0]*(weak.DR@weak.shape(entropy)@weak.Z.T)+j[...,2]*(weak.R@weak.shape(entropy)@weak.DZ.T))
                    lumped=adv/d.volume.reshape(-1,d.nz)
                    exact=np.linalg.solve(mr,np.linalg.solve(mz,adv.T).T)
                    baseline={k:np.array_equal(value,previous[label+'_'+k]) for k,value in [('original',original),('lumped',lumped),('exact_mass',exact)]}
                    baseline_comparisons={k:relative(value,previous[label+'_'+k]) for k,value in [('original',original),('lumped',lumped),('exact_mass',exact)]} if elements==8 else None
                    if elements==8 and not all(baseline.values()):raise ValueError('All original F8 coefficient vectors must match exactly')
                    item=dict(name=label,baseline_coefficients_exact=baseline if elements==8 else None,
                        baseline_coefficient_comparisons=baseline_comparisons,
                        reconstructed_source_vs_fixed=compare(sample,reconstructed_source,fixed_source),
                        errors={},conservation={})
                    for kind,value in [('original',original),('lumped',lumped),('exact_mass',exact)]:
                        physical=sample.scalar(value.reshape(d.shape[:-1]))/n
                        item['errors'][kind]=compare(sample,physical,target)
                        item['errors'][kind]['squared_integral']=float(np.sum(sample.weight*physical**2))
                        item['errors'][kind]['fingerprint']=fingerprint(physical)
                        terms=sample.weight*physical*n
                        item['conservation'][kind]=dict(physical_integral=float(terms.sum()),absolute_integral=float(abs(terms).sum()))
                        arrays[label+'_'+kind]=value
                    row['species'].append(item)
                path=destination/f'F{elements}_{name}_q{order}.npz';np.savez_compressed(path,**arrays)
                row.update(arrays_file=path.name,arrays_sha256=sha(path));result['cases'].append(row);save()
                print('fixed-current entropy space',order,elements,name,flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.species_space_audit'])!=identity:
        raise ValueError('Entropy-space inputs/source changed')
    result.update(terminal=True,all_cases_complete=True,
                  scope='Only the entropy source representation is refined; no complete coupled grid-convergence, growth-causality or formation claim.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
