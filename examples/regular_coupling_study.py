"""Original force/current checks on the shared regular-coordinate field map."""
import argparse
import json
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from plasma_bench.variable_compatible_current import weak_vector
from plasma_bench.kernel_induction import KernelInduction
from .regular_coupling_reference import initialize,force,force_regular


def relative(diff,target,weight):
    if diff.ndim==weight.ndim+1:weight=weight[...,None]
    return float(np.sqrt(np.sum(weight*diff*diff)/np.sum(weight*target*target)))


def run_case(spec):
    start=time.monotonic();ref,d,u,m,field=initialize(spec);g=m.geometry;mapping=m.mapping
    physical_force=force(g.r,g.z,ref);forced=mapping.force_rate(u,physical_force)
    exact=force_regular(d.r,d.z,ref)/d.mass;exact[0,0,:,:2]=0
    error=forced[...,1:4]-exact
    f=d.fields(u);work=float(np.sum(d.volume[...,None]*f['energy_gradient']*forced));physical_work=float(np.sum(g.weight*np.sum(mapping.velocity(u)*physical_force,axis=-1)))
    axial=float(d.mass*np.sum(d.volume*forced[...,3]));torque=float(d.mass*np.sum(d.volume*2*d.x*forced[...,2]))
    axial_exact=float(np.sum(g.weight*physical_force[...,2]));torque_exact=float(np.sum(g.weight*g.r*physical_force[...,1]))
    full_rate,_=d.derivative(u);nt=full_rate[...,0];nt_sample=m.density_map.sample(nt);analytic=ref.fields(g.r,g.z)
    current=field.current(m);j=g.sample_electric(current)
    other=weak_vector(g,analytic['other_forcing']);density_force=m.density_forcing(current,nt)
    solved=KernelInduction(m,tolerance=1e-10,max_iterations=3000).solve(other+density_force)
    psi,gdot=g.unpack(solved.magnetic);bt=g.magnetic_samples(psi,gdot);jt=g.sample_electric(solved.current)
    kr=m.current_energy_rate(current,solved.current,nt);br=m.magnetic_energy_rate(field,solved.magnetic)
    total=float(np.sum(m.mesh.volumes*(kr['total']+br)))
    target=float(-current@other-np.sum(m.mesh.volumes*kr['density']))
    scale=float(np.sum(g.weight*(abs(np.sum(field.samples(g)*bt,axis=-1)/MU0)+abs(m.inertia_coefficient*np.sum(j*jt,axis=-1))+abs(.5*m.inertia_coefficient*np.sum(j*j,axis=-1)*nt_sample/m.density))))
    # Independent complex differentiation of the physical integral at fixed J.
    dn=u[...,0]*(.2+.1*np.cos(d.z*np.pi/ref.half_length));dnh=m.density_map.sample(dn)
    h=1e-24;complex_energy=np.sum(g.weight*.5*m.reduced_mass*np.sum(j*j,axis=-1)/(E**2*(m.density+1j*h*dnh)))
    gradient_work=float(np.sum(m.density_energy_gradient(current)*dn));complex_work=float(np.imag(complex_energy)/h)
    pair=m.electromagnetic_pair(u,field,nt);coupled=KernelInduction(m,tolerance=1e-10,max_iterations=3000).solve(pair['forcing'])
    node_fields=ref.fields(d.r,d.z);node_b=node_fields['magnetic'].copy();node_b[...,2]+=ref.axial_bias
    node_force=np.cross(node_fields['current'],node_b)/d.mass;lorentz_target=node_force.copy()
    np.divide(node_force[...,:2],d.r[...,None],out=lorentz_target[...,:2],where=d.r[...,None]>0);lorentz_target[0,0,:,:2]=0
    cr=m.current_energy_rate(current,coupled.current,nt);cb=m.magnetic_energy_rate(field,coupled.magnetic)
    bulk=float(np.sum(d.volume[...,None]*f['energy_gradient']*pair['fluid_rate']))
    coupled_total=float(bulk+np.sum(m.mesh.volumes*(cr['total']+cb)))
    coupled_scale=abs(bulk)+float(np.sum(m.mesh.volumes*(abs(cr['total'])+abs(cb))))+abs(pair['unresolved_density_exchange'])
    record=dict(spec=spec,status='complete',force_relative_L2=[relative(error[...,k],exact[...,k],d.volume) for k in range(3)],
        axis_force_relative_L2=[relative(error[0,1:,:,k],exact[0,1:,:,k],d.volume[0,1:]) for k in range(3)],
        force_work_backward=abs(work-physical_work)/max(abs(work),abs(physical_work)),
        axial_force_relative=abs(axial-axial_exact)/abs(axial_exact),angular_torque_relative=abs(torque-torque_exact)/abs(torque_exact),
        current_error=relative(j-analytic['current'],analytic['current'],g.weight),
        magnetic_rate_error=relative(bt-analytic['magnetic_rate'],analytic['magnetic_rate'],g.weight),
        current_rate_error=relative(jt-analytic['current_rate'],analytic['current_rate'],g.weight),
        continuity_error=relative(nt_sample-analytic['nt'],analytic['nt'],g.weight),
        lorentz_force_error=[relative(pair['fluid_rate'][...,k+1]-lorentz_target[...,k],lorentz_target[...,k],d.volume) for k in range(3)],
        axis_lorentz_force_error=[relative(pair['fluid_rate'][0,1:,:,k+1]-lorentz_target[0,1:,:,k],lorentz_target[0,1:,:,k],d.volume[0,1:]) for k in range(3)],
        induction_energy_backward=abs(total-target)/scale,density_gradient_relative=abs(gradient_work-complex_work)/abs(complex_work),
        ideal_work_backward=abs(bulk-current@pair['ideal_forcing'])/max(abs(bulk),abs(current@pair['ideal_forcing'])),
        hall_work_backward=float(abs(current@pair['hall_forcing'])/(abs(current)@abs(pair['hall_forcing']))),
        incomplete_total_energy_rate=coupled_total,unresolved_density_exchange=pair['unresolved_density_exchange'],
        exchange_balance_backward=abs(coupled_total-pair['unresolved_density_exchange'])/coupled_scale,
        residual=solved.residual,residual_components=solved.residual_components,coupled_residual=coupled.residual,
        elapsed_seconds=time.monotonic()-start)
    arrays=dict(state=u,fluid_rate=full_rate,force_rate=forced,exact_force_rate=exact,magnetic=field.packed_increment(g),current=current,
        magnetic_rate=solved.magnetic,current_rate=solved.current,electric=solved.electric,
        density_gradient=m.density_energy_gradient(current),density_direction=dn,coupled_magnetic_rate=coupled.magnetic,
        coupled_current_rate=coupled.current,coupled_electric=coupled.electric,lorentz_rate=pair['fluid_rate'],
        ideal_forcing=pair['ideal_forcing'],hall_forcing=pair['hall_forcing'],density_forcing=pair['density_forcing'])
    return record,arrays


def protocol():
    specs=[dict(elements=ne,degree=p) for p in (3,5) for ne in (2,4,8,16)]
    specs.extend([dict(elements=8,quadrature=q) for q in (8,16)])
    specs.extend([dict(elements=8,contrast=0.),dict(elements=8,density=1e16),dict(elements=8,density=1e20)])
    specs.extend([dict(elements=8,field_nr=8),dict(elements=8,field_nr=32),dict(elements=8,points=80)])
    return specs


def main():
    from .regular_coupling_records import provenance,freeze,sha
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--cases',type=int,nargs='*');args=p.parse_args()
    directory=args.directory;identity=provenance()
    if directory.exists():
        if json.loads((directory/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source mismatch')
    else:directory.mkdir(parents=True);freeze(directory,identity)
    specs=protocol()
    for i in args.cases if args.cases is not None else range(len(specs)):
        path=directory/f'case_{i:02d}.json'
        if path.exists():raise ValueError('Terminal evidence already exists')
        try:
            result,arrays=run_case(specs[i]);raw=path.with_suffix('.npz');np.savez_compressed(raw,**arrays)
            result.update(arrays_file=raw.name,arrays_sha256=sha(raw))
        except Exception as exc:result=dict(spec=specs[i],status='failed',error=repr(exc),traceback=traceback.format_exc())
        result['source_sha256']=identity['source_sha256'];path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
        if result['status']!='complete':raise RuntimeError('Failure preserved; no retry')
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')


if __name__=='__main__':main()
