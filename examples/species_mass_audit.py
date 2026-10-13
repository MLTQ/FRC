"""Independent mass-map and represented-source entropy diagnostics."""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.constants import elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .species_weak_sampling import WeakSampler
from .species_interface_audit import relative
from .trajectory_tangent_probe import verified_inputs
from .trajectory_tangent_audit import compare
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def mass_reference(fluid):
    vr=fluid.volume.reshape(-1,fluid.nz)[:,0]/fluid.dz
    radial=np.zeros((len(vr),fluid.elements))
    for e in range(fluid.elements):
        nodes=fluid.x[e,:,0];xi=2*(nodes-nodes[0])/(nodes[-1]-nodes[0])-1
        rows=slice(e*(fluid.degree+1),(e+1)*(fluid.degree+1))
        mode=np.sqrt(vr[rows])*Legendre.basis(fluid.degree)(xi)
        radial[rows,e]=mode/np.linalg.norm(mode)
    axial=(-1.)**np.arange(fluid.nz);axial/=np.linalg.norm(axial)
    radial_ratio=fluid.degree/(2*fluid.degree+1)
    return (np.eye(len(vr))-(1-radial_ratio)*(radial@radial.T),
            np.eye(fluid.nz)-.5*np.outer(axial,axial),radial,axial,vr)


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh mass-audit output required')
    prior=Path('outputs/species_interface_audit')
    record,inputs=verified_inputs(prior,str(prior)+'_replay_verification')
    inputs.update(record['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original replayed weak audit required')
    with np.load('outputs/trajectory_tangent_probe/context.npz') as raw:context=dict(raw)
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    d=RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
                             degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz']);base=context['state']
    mr_reference,mz_reference,pr,pz,vr=mass_reference(d)
    alpha=ref.electron_mass/(ref.ion_mass+ref.electron_mass);beta=ref.ion_mass/(ref.ion_mass+ref.electron_mass)
    identity=provenance(['examples.species_mass_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,all_cases_complete=False,source_sha256=identity['source_sha256'],input_sha256=inputs,
                spec=spec,quadrature=[],analytic_radial_mass_ratio=5/11,analytic_axial_nyquist_mass_ratio=.5)
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save()
    for order in (16,20):
        weak=WeakSampler(d,mesh,order);sample=weak.sample;n=sample.scalar(base[...,0])
        mr=weak.R.T@(weak.wr[:,None]*weak.R);mz=weak.Z.T@(weak.wz[:,None]*weak.Z)
        nmr=mr/np.sqrt(vr[:,None]*vr[None,:]);nmz=mz/d.dz
        q=dict(order=order,cases=[],radial_mass_vs_analytic=relative(nmr,mr_reference),
               axial_mass_vs_analytic=relative(nmz,mz_reference),
               radial_eigenvalues=np.linalg.eigvalsh(nmr).tolist(),axial_eigenvalues=np.linalg.eigvalsh(nmz).tolist())
        result['quadrature'].append(q)
        for name in ('saved_error','independent_mode'):
            j=sample.electric(context[name+'_current'])
            with np.load(prior/f'{name}_q{order}.npz') as raw:old=dict(raw)
            saved=dict(radial_mass=mr,axial_mass=mz,radial_analytic=mr_reference,axial_analytic=mz_reference)
            row=dict(direction=name,species=[])
            for index,(label,a) in enumerate(zip(('electron','ion'),(-beta,alpha))):
                slot=index+4;s0=base[...,slot]/base[...,0]
                source=-a/E*(j[...,0]*sample.scalar(s0,axis=0)+j[...,2]*sample.scalar(s0,axis=2))
                target=source/n
                adv=-a/E*old[label+'_production_advection']
                original=old[label+'_production_Cdot'].reshape(-1,d.nz)
                lumped=adv/d.volume.reshape(-1,d.nz)
                exact=np.linalg.solve(mr,np.linalg.solve(mz,adv.T).T)
                independently_integrated=weak.adjoint(source)
                reconstructed={k:sample.scalar(v.reshape(d.shape[:-1]))/n for k,v in
                               [('original',original),('lumped',lumped),('exact_mass',exact)]}
                components=dict(interface_transfer=reconstructed['original']-reconstructed['lumped'],
                                mass_map=reconstructed['lumped']-reconstructed['exact_mass'],
                                representation=reconstructed['exact_mass']-target)
                total=reconstructed['original']-target
                # Store all actual signed products; cancellation cannot be inferred from norms.
                names=list(components);products=np.array([[np.sum(sample.weight*components[x]*components[y]) for y in names] for x in names])
                target_norm=float(np.sum(sample.weight*target**2));difference_norm=float(np.sum(sample.weight*total**2))
                normcoords=np.sqrt(vr[:,None]*d.dz)*exact
                highr=pr@(pr.T@normcoords);highz=(normcoords@pz)[:,None]*pz[None,:]
                highboth=pr@(pr.T@highz)
                modes={'lower_radial_lower_axial':normcoords-highr-highz+highboth,
                       'highest_radial_lower_axial':highr-highboth,
                       'lower_radial_nyquist':highz-highboth,'highest_radial_nyquist':highboth}
                mode_squares={k:float(np.sum(v*v)) for k,v in modes.items()};mode_total=float(np.sum(normcoords**2))
                conservation={}
                for kind,value in [('original',original),('lumped',lumped),('exact_mass',exact)]:
                    native=d.volume.reshape(-1,d.nz)*value
                    physical=sample.weight*sample.scalar(value.reshape(d.shape[:-1]))
                    conservation[kind]=dict(native_integral=float(native.sum()),native_absolute_integral=float(abs(native).sum()),
                                             physical_integral=float(physical.sum()),physical_absolute_integral=float(abs(physical).sum()))
                item=dict(name=label,advective_covector_vs_fresh_integral=relative(adv,independently_integrated),
                    field_identity=compare(sample,sum(components.values()),total),
                    strong_source_squared=target_norm,actual_difference_squared=difference_norm,
                    rate_comparisons={k:compare(sample,v,target) for k,v in reconstructed.items()},
                    component_names=names,component_products=products.tolist(),
                    component_norms_relative_to_strong={k:float(np.sqrt(np.sum(sample.weight*v*v)/target_norm)) for k,v in components.items()},
                    components_fingerprints={k:fingerprint(v) for k,v in components.items()},
                    squared_norm_sum_defect=float(abs(products.sum()-difference_norm)/difference_norm),
                    point_volume_normalized_mode_fractions={k:v/mode_total for k,v in mode_squares.items()},
                    mode_partition_relative_defect=float(abs(sum(mode_squares.values())-mode_total)/mode_total),
                    conservation=conservation)
                row['species'].append(item)
                for kind,value in [('original',original),('lumped',lumped),('exact_mass',exact),('advective_covector',adv)]:
                    saved[label+'_'+kind]=value
                saved[label+'_products']=products
                saved[label+'_mode_squares']=np.array(list(mode_squares.values()))
            path=destination/f'{name}_q{order}.npz';np.savez_compressed(path,**saved)
            row.update(arrays_file=path.name,arrays_sha256=sha(path));q['cases'].append(row);save()
            print('physical entropy mass budget',order,name,flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.species_mass_audit'])!=identity:
        raise ValueError('Original mass-audit source/input required')
    result.update(terminal=True,all_cases_complete=True,
                  scope='Diagnostic exact physical mass and representation split only; no evolved rate or coupled-growth causality claim.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
