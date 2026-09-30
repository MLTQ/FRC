"""Measured original current transport and drift stress with face controls."""
import argparse
import json
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.current_transport import CurrentTransport
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.physical_derivatives import TensorSurface,CurrentDifferential,fluid_map
from .regular_coupling_reference import initialize
from .regular_coupling_study import relative
from .current_transport_reference import exact


def run_case(spec):
    start=time.monotonic();ref,d,u,m,field=initialize(spec);g=m.geometry;current=field.current(m)
    fluid_rate,_=d.derivative(u);nt=fluid_rate[...,0];transport=CurrentTransport(m);pair=transport.pair(u,current,nt)
    oracle=exact(ref,g.r,g.z,m.delta);electric_errors={};electric_axis={};electric_max={};values={}
    for key in ('advection','stretching','self_drift'):
        coefficients=np.concatenate([solve(a) for solve,a in zip(g.solves,g.split_electric(pair[key]))]);value=g.sample_electric(coefficients);values[key]=coefficients
        wanted=oracle[key];electric_errors[key]=relative(value-wanted,wanted,g.weight)
        electric_axis[key]=relative((value-wanted)[0],wanted[0],g.weight[0])
        electric_max[key]=float(np.max(np.linalg.norm(value-wanted,axis=-1))/np.sqrt(np.sum(g.weight[...,None]*wanted*wanted)/np.sum(g.weight)))
    target=exact(ref,d.r,d.z,m.delta)['force']/d.mass
    np.divide(target[...,:2],d.r[...,None],out=target[...,:2],where=d.r[...,None]>0);target[0,0,:,:2]=0
    error=pair['fluid_rate'][...,1:4]-target
    j=g.sample_electric(current);v=m.mapping.velocity(u);n=m.density;vg=transport.fluid.velocity_gradient(u);jg=transport.current.gradient(current)
    stress=m.reduced_mass/E**2*j[..., :,None]*j[...,None,:]/n[...,None,None]
    omitted_force=transport.fluid.stress_rate(u,stress)
    omitted_adv=transport.advection(current,v,nt,[]);omitted_coeff=np.concatenate([solve(a) for solve,a in zip(g.solves,g.split_electric(omitted_adv))])
    nt_sample=m.density_map.sample(nt);continuity=nt_sample+n*np.trace(vg,axis1=-2,axis2=-1)+v[...,0]*transport.fluid.scalar(u[...,0],0)+v[...,2]*transport.fluid.scalar(u[...,0],2)
    divj=np.trace(jg,axis1=-2,axis2=-1);jnorm=np.sum(g.weight*np.sum(j*j,axis=-1));jumps=[]
    for radial,surface,maps,fluids in transport.faces:
        normal=0 if radial else 2;jl,jr=(cm.sample(current) for cm in maps)
        jumps.append(float(np.sqrt(ref.radius*np.sum(surface.weight*(jl[...,normal]-jr[...,normal])**2)/jnorm)))
    wall=TensorSurface(g,True);wall.radial_points=np.array([ref.radius]);wall.r,wall.z=np.meshgrid(wall.radial_points,wall.axial_points,indexing='ij');wall.weight=2*np.pi*wall.r*g.wz[None,:]
    wall_j=CurrentDifferential(wall,radial_cells=np.array([m.mesh.nr-1])).sample(current)
    wall_v=fluid_map(d,wall,'left').velocity(u)
    wall_normal=float(np.sqrt(ref.radius*np.sum(wall.weight*wall_j[...,0]**2)/jnorm))
    gradient=d.fields(u)['energy_gradient'];bulk=float(np.sum(d.volume[...,None]*gradient*pair['fluid_rate']))
    density=float(-np.sum(m.density_energy_gradient(current)*nt));adv=float(current@pair['advection']);stretch=float(current@pair['stretching']);self_work=float(current@pair['self_drift'])
    base=m.electromagnetic_pair(u,field,nt);solved=KernelInduction(m,tolerance=1e-10,max_iterations=3000).solve(base['forcing']+pair['forcing'])
    all_fluid=fluid_rate+base['fluid_rate']+pair['fluid_rate'];kinetic_thermal=float(np.sum(d.volume[...,None]*gradient*all_fluid))
    kr=m.current_energy_rate(current,solved.current,nt);br=m.magnetic_energy_rate(field,solved.magnetic)
    total=kinetic_thermal+float(np.sum(m.mesh.volumes*(kr['total']+br)))
    scale=float(np.sum(d.volume[...,None]*abs(gradient*all_fluid))+np.sum(m.mesh.volumes*(abs(kr['total'])+abs(br))))
    backward=lambda a,b:abs(a-b)/max(abs(a),abs(b),np.finfo(float).tiny)
    record=dict(spec=spec,status='complete',electric_relative_L2=electric_errors,electric_first_field_cell=electric_axis,electric_maximum_over_rms=electric_max,
        force_relative_L2=[relative(error[...,k],target[...,k],d.volume) for k in range(3)],
        force_first_fluid_element=[relative(error[0,1:,:,k],target[0,1:,:,k],d.volume[0,1:]) for k in range(3)],
        force_without_faces=[relative(omitted_force[...,k+1]-target[...,k],target[...,k],d.volume) for k in range(3)],
        advection_without_faces=relative(g.sample_electric(omitted_coeff)-oracle['advection'],oracle['advection'],g.weight),
        continuity_point_defect=relative(continuity,ref.fields(g.r,g.z)['nt'],g.weight),
        charge_divergence_scaled=float(np.sqrt(ref.radius**2*np.sum(g.weight*divj*divj)/jnorm)),normal_current_jumps=jumps,
        wall_normal_current_scaled=wall_normal,wall_radial_velocity_max=float(np.max(abs(wall_v[...,0]))),
        advection_work=adv,density_exchange=density,advection_work_backward=backward(adv,density),stress_work=bulk,stretching_work=stretch,
        stress_work_backward=backward(bulk,stretch),self_drift_work=self_work,
        self_drift_work_backward=float(abs(self_work)/(abs(current)@abs(pair['self_drift']))),
        axial_force_scaled=float(abs(d.mass*np.sum(d.volume*pair['fluid_rate'][...,3]))/(d.mass*np.sum(d.volume*abs(pair['fluid_rate'][...,3])))),
        angular_torque_scaled=float(abs(d.mass*np.sum(d.volume*2*d.x*pair['fluid_rate'][...,2]))/(d.mass*np.sum(d.volume*2*d.x*abs(pair['fluid_rate'][...,2])))),
        total_energy_rate=total,total_energy_backward=abs(total)/scale,residual=solved.residual,residual_components=solved.residual_components,elapsed_seconds=time.monotonic()-start)
    arrays=dict(state=u,fluid_rate=fluid_rate,full_fluid_rate=all_fluid,magnetic=field.packed_increment(g),current=current,
        magnetic_rate=solved.magnetic,current_rate=solved.current,electric=solved.electric,stress_rate=pair['fluid_rate'],exact_stress_rate=target,
        advection=pair['advection'],stretching=pair['stretching'],self_drift=pair['self_drift'],omitted_stress_rate=omitted_force,omitted_advection=omitted_adv,
        **{key+'_electric':value for key,value in values.items()})
    return record,arrays


def protocol():
    specs=[dict(elements=n,degree=p) for p in (3,5) for n in (4,8,16,32)]
    specs.extend([dict(elements=8,quadrature=q) for q in (8,16)])
    specs.extend([dict(elements=8,field_nr=n) for n in (8,32)])
    specs.extend([dict(elements=8,points=80),dict(elements=8,contrast=0.),dict(elements=8,density=1e16),dict(elements=8,density=1e20)])
    return specs


def main():
    from .current_transport_records import provenance,freeze,sha
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--cases',type=int,nargs='*');args=p.parse_args();directory=args.directory;identity=provenance()
    if directory.exists():
        if json.loads((directory/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source mismatch')
    else:directory.mkdir(parents=True);freeze(directory,identity)
    specs=protocol()
    for i in args.cases if args.cases is not None else range(len(specs)):
        path=directory/f'case_{i:02d}.json'
        if path.exists():raise ValueError('Terminal evidence already exists')
        try:
            r,a=run_case(specs[i]);raw=path.with_suffix('.npz');np.savez_compressed(raw,**a);r.update(arrays_file=raw.name,arrays_sha256=sha(raw))
        except Exception as exc:r=dict(spec=specs[i],status='failed',error=repr(exc),traceback=traceback.format_exc())
        r['source_sha256']=identity['source_sha256'];path.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
        if r['status']!='complete':raise RuntimeError('Failure preserved; no retry')
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')


if __name__=='__main__':main()
