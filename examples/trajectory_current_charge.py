"""Independent one-sided current traces and cellwise cylindrical Gauss audit."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import numpy as np
from numpy.polynomial import Polynomial
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler, bases
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def face_edges(edges, periodic=False):
    count=len(edges)-1;result=np.zeros((2*count,count))
    for cell in range(count):
        start=cell-2
        if not periodic:start=min(max(start,0),count-5)
        ids=np.arange(start,start+6)
        nodes=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
        scale=nodes[-1]-nodes[0]
        for side in (0,1):
            roots=(nodes-edges[cell+side])/scale;derivatives=[]
            for i in range(6):
                other=np.delete(roots,i);p=Polynomial.fromroots(other)/np.prod(roots[i]-other)
                derivatives.append(p.deriv()(0.)/scale)
            for i,index in enumerate(ids[:-1]):
                result[2*cell+side,index%count]=-sum(derivatives[:i+1])*(nodes[i+1]-nodes[i])
    return result


def integrate_faces(mesh, sample, radial, axial, weights):
    order=len(weights);rface=np.column_stack((mesh.rf[:-1],mesh.rf[1:]))
    rweight=2*np.pi*rface[:,:,None,None]*mesh.dz*.5*weights[None,None,None,:]
    zweight=np.pi*mesh.dr*sample.r.reshape(mesh.nr,order)[:,:,None,None]*weights[None,:,None,None]
    rf=np.sum(radial*rweight,axis=-1)
    zf=np.sum(axial*zweight,axis=1)
    net=rf[:,1]-rf[:,0]+zf[:,:,1]-zf[:,:,0]
    bound=abs(rf[:,1])+abs(rf[:,0])+abs(zf[:,:,1])+abs(zf[:,:,0])
    return rf,zf,net,bound,rweight,zweight


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh charge-audit output required')
    if shutil.disk_usage('.').free<384*1024**2:raise RuntimeError('128MiB payload plus256MiB free reserve required')
    prior=Path('outputs/trajectory_spatial_parts_audit')
    verifier=Path('outputs/trajectory_spatial_parts_audit_replay_verification/summary.json')
    verification=json.loads(verifier.read_text());audit=json.loads((prior/'summary.json').read_text())
    if not verification['terminal'] or not verification['bitwise_metrics'] or sha(prior/'summary.json')!=verification['original_summary_sha256']:
        raise ValueError('Complete original physical contribution replay required')
    inputs={str(verifier):sha(verifier),str(prior/'summary.json'):sha(prior/'summary.json')}
    prior_source=json.loads((prior/'provenance.json').read_text())
    if any(sha(p)!=h for p,h in prior_source['source_files'].items()):raise ValueError('Original physical audit source required')
    inputs[str(prior/'provenance.json')]=sha(prior/'provenance.json')
    inputs.update(audit['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Unchanged complete audit inputs required')
    probe=Path('outputs/trajectory_tangent_probe');record=json.loads((probe/'summary.json').read_text())
    source=json.loads((probe/'provenance.json').read_text())
    if any(sha(p)!=h for p,h in source['source_files'].items()):raise ValueError('Original tangent source required')
    with np.load(probe/'context.npz') as raw:context=dict(raw)
    mode=ContinuumContactMode();ref=mode.ref;spec=record['spec']
    d=RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
                            degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    identity=provenance(['examples.trajectory_current_charge']);destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(terminal=False,source_sha256=identity['source_sha256'],input_sha256=inputs,spec=spec,cases=[])
    resources=[]
    def observe():
        resources.append(dict(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,free_bytes=shutil.disk_usage('.').free))
        if resources[-1]['peak_rss_bytes']>8*1024**3 or resources[-1]['free_bytes']<256*1024**2:
            raise RuntimeError('Recorded charge-audit resource bounds exceeded')
    def save():
        p=destination/'summary.pending.json';p.write_text(json.dumps(result,indent=2)+'\n');p.replace(destination/'summary.json')
    save();observe()
    ex=face_edges(.5*mesh.rf**2);ez=face_edges(mesh.zf,True)
    rf=np.column_stack((mesh.rf[:-1],mesh.rf[1:])).ravel()
    zf=np.column_stack((mesh.zf[:-1],mesh.zf[1:])).ravel()
    count=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    for order in (16,20):
        sample=TransportSampler(d,mesh,order);_,w=np.polynomial.legendre.leggauss(order)
        nx,_,_,_,_=bases(.5*mesh.rf**2,.5*sample.r**2,5)
        nz,_,_,_,_=bases(mesh.zf,sample.z,5,True)
        # Separate analytic quadrature control, with div J exactly zero.
        k=2*np.pi/ref.length
        ar=(rf[:,None]*k*np.cos(k*sample.z[None,:])).reshape(mesh.nr,2,mesh.nz,order)
        az=np.broadcast_to(-2*np.sin(k*zf)[None,:],(len(sample.r),len(zf))).reshape(mesh.nr,order,mesh.nz,2)
        _,_,analytic_net,analytic_bound,_,_=integrate_faces(mesh,sample,ar,az,w)
        control=dict(max_cell_flux_A=float(np.max(abs(analytic_net))),
                     zero_divergence_flux_relative=float(np.linalg.norm(analytic_net)/np.linalg.norm(analytic_bound)))
        for name in ('saved_error','independent_mode'):
            current=context[name+'_current'];j,gradient=sample.electric(current,gradient=True)
            div=np.trace(gradient,axis1=-2,axis2=-1)
            volume=np.sum((sample.weight*div).reshape(mesh.nr,order,mesh.nz,order),axis=(1,3))
            jr=current[:count].reshape(mesh.nr,mesh.nz);jz=current[count+theta:].reshape(mesh.nr+1,mesh.nz)
            radial=(rf[:,None]*(ex@jr@nz.T)).reshape(mesh.nr,2,mesh.nz,order)
            axial=(nx@jz@ez.T).reshape(mesh.nr,order,mesh.nz,2)
            rflux,zflux,net,bound,rweight,zweight=integrate_faces(mesh,sample,radial,axial,w)
            jump_r=radial[:-1,1]-radial[1:,0]
            jump_z=axial[:,:,:,1]-np.roll(axial[:,:,:,0],-1,axis=2)
            ri=float(np.sum(jump_r**2*rweight[:-1,1]));zi=float(np.sum(jump_z**2*zweight[...,0]))
            wall=float(np.sum(radial[-1,1]**2*rweight[-1,1]))
            ji=float(np.sum(sample.weight*np.sum(j*j,axis=-1)));di=float(np.sum(sample.weight*(div*div)))
            previous=next(c for q in audit['quadrature'] if q['order']==order for c in q['cases'] if c['direction']==name)
            if ji!=previous['interior_diagnostics']['current_squared_integral'] or di!=previous['interior_diagnostics']['current_divergence_squared_integral']:
                raise ValueError('Exact previous physical divergence integrals required')
            norm_bound=float(np.linalg.norm(bound));defect=float(np.linalg.norm(volume-net))
            row=dict(order=order,direction=name,current_squared_integral=ji,divergence_squared_integral=di,
                radius_times_divergence_over_current=float(ref.radius*np.sqrt(di/ji)),
                radial_jump_squared_surface_integral=ri,axial_jump_squared_surface_integral=zi,
                wall_normal_squared_surface_integral=wall,
                normalized_radial_jump=float(np.sqrt(ref.radius*ri/ji)),normalized_axial_jump=float(np.sqrt(ref.radius*zi/ji)),
                normalized_wall_normal=float(np.sqrt(ref.radius*wall/ji)),
                max_radial_jump_A_m2=float(np.max(abs(jump_r))),max_axial_jump_A_m2=float(np.max(abs(jump_z))),
                gauss_relative_to_absolute_face_flux=defect/norm_bound,
                gauss_max_cell_defect_A=float(np.max(abs(volume-net))),
                summed_volume_divergence_A=float(volume.sum()),summed_cell_outward_flux_A=float(net.sum()),
                wall_outward_flux_A=float(rflux[-1,1].sum()),axis_outward_flux_A=float(-rflux[0,0].sum()),
                analytic_solenoidal_control=control,current_fingerprint=fingerprint(j),divergence_fingerprint=fingerprint(div))
            path=destination/f'{name}_q{order}.npz'
            np.savez_compressed(path,current_coefficients=current,r=sample.r,z=sample.z,divergence=div,
                radial_traces=radial,axial_traces=axial,radial_jumps=jump_r,axial_jumps=jump_z,
                cell_volume_divergence=volume,cell_radial_flux=rflux,cell_axial_flux=zflux,
                cell_outward_flux=net,cell_absolute_face_flux_bound=bound,analytic_cell_flux=analytic_net,
                analytic_cell_bound=analytic_bound)
            row.update(arrays_file=path.name,arrays_sha256=sha(path));result['cases'].append(row);save();observe()
            print('one-sided current charge',name,order,row['normalized_radial_jump'],row['normalized_axial_jump'],flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.trajectory_current_charge'])!=identity:
        raise ValueError('Charge-audit input or source changed')
    resource_path=destination/'resources.json';resource_path.write_text(json.dumps(dict(stage_observations=resources,
        continuous_guard=False,maximum_rss_bytes=8*1024**3,minimum_free_bytes=256*1024**2),indent=2)+'\n')
    result.update(terminal=True,all_cases_complete=True,
        scope='Volume and one-sided physical charge-flux consistency. No replacement current, discrete-graph identification or causal instability claim. Normalized directions, not physical trajectory amplitudes.')
    save();return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
