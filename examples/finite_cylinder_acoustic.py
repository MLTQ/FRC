"""Complete finite-tube acoustic evolution with physical end-pressure forces."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, electron_mass as ME, elementary_charge as E
from plasma_bench.finite_cylinder_fluid import FiniteCylinderFluid
from plasma_bench.regular_field_map import polynomial
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, *, axial_points=12, steps=1600, amplitude=1e-5, mode=1):
    destination=Path(destination)
    if destination.exists() or axial_points not in (8,12,16,24) or type(steps) is not int or steps<4:
        raise ValueError('Fresh output,8/12/16/24 axial points and positive full-period step count required')
    if not np.isfinite(amplitude) or not 0<amplitude<.01:raise ValueError('Explicit small positive acoustic amplitude required')
    if type(mode) is not int or mode<1:raise ValueError('Positive integer acoustic mode required')
    started=time.monotonic();identity=provenance(['examples.finite_cylinder_acoustic'])
    destination.mkdir(parents=True);freeze(destination,identity)
    radius=.075;length=1.;density=1e20;mass=39.948*atomic_mass+ME;temperature=np.array([5.,1.]);gamma=5/3
    fluid=FiniteCylinderFluid(radius,length,2,axial_points,mass,degree=5,gammas=(gamma,gamma))
    pressure=density*E*temperature.sum();speed=np.sqrt(gamma*pressure/(mass*density));wave=mode*np.pi/length
    duration=2*length/(mode*speed);dt=duration/steps
    entropy=np.log(temperature)-(gamma-1)*np.log(density/fluid.thermal.n0)
    state=np.zeros(fluid.shape);state[...,0]=density*(1+amplitude*np.cos(wave*(fluid.z+.5*length)))
    state[...,4:]=state[...,0,None]*entropy
    initial_state=state.copy();initial_fields=fluid.fields(state);volume=np.pi*radius**2*length
    # Independent denser axial root-product reconstruction; actual polynomial
    # field, not a periodic interpolant. Radial target is exactly uniform.
    q,w=np.polynomial.legendre.leggauss(64);z=.5*length*q;weights=.5*w
    mapping=polynomial(fluid.z.ravel(),z)
    radial=fluid.radial_weight[...,0]/np.sum(fluid.radial_weight[...,0])
    sample=lambda value:np.einsum('er,erj,ij->i',radial,value,mapping)
    norm=lambda value:float(np.sqrt(np.sum(weights*value*value)))
    initial_n=sample(initial_state[...,0])
    initial_s=np.stack([sample(initial_state[...,k])/initial_n for k in (4,5)],axis=-1)
    initial_t=np.exp(initial_s+(gamma-1)*np.log(initial_n[:,None]/fluid.thermal.n0))
    def stable_internal_change(n,entropy_density,base_n,base_s,base_h):
        dn=n-base_n;ratio=dn/base_n
        ds=(entropy_density-base_n[...,None]*base_s-base_s*dn[...,None])/n[...,None]
        log_change=ds+(gamma-1)*np.log1p(ratio)[...,None]
        return np.sum(base_n[...,None]*base_h*(ratio[...,None]+(1+ratio[...,None])*np.expm1(log_change)),axis=-1)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,
        radius_m=radius,length_m=length,number_density_m3=density,mass_kg=mass,initial_temperature_ev=temperature.tolist(),
        gamma=gamma,acoustic_mode=mode,axial_points=axial_points,radial_elements=2,radial_degree=5,amplitude=amplitude,
        sound_speed_m_s=speed,source_sha256=identity['source_sha256'],records=[],
        scope='Full nonlinear two-entropy Euler in a sealed finite cylinder. One complete acoustic period against the independent linear reflected-wave limit. No magnetic/current coupling or formation claim.')

    def observe(u,when):
        f=fluid.fields(u);rate,diagnostic=fluid.derivative(u)
        n=sample(u[...,0]);mz=sample(u[...,3]);v=mz/n
        specific=np.stack([sample(u[...,k])/n for k in (4,5)],axis=-1)
        t=np.exp(specific+(gamma-1)*np.log(n[:,None]/fluid.thermal.n0));p=n*E*t.sum(axis=1)
        phase=wave*(z+.5*length);omega=wave*speed
        expected_n=density*(1+amplitude*np.cos(phase)*np.cos(omega*when))
        expected_v=amplitude*speed*np.sin(phase)*np.sin(omega*when)
        expected_p=pressure*(1+gamma*amplitude*np.cos(phase)*np.cos(omega*when))
        kinetic=.5*mass*n*v*v;internal=n*E*t.sum(axis=1)/(gamma-1)
        native_change=stable_internal_change(u[...,0],u[...,4:],initial_state[...,0],initial_fields['entropy'],initial_fields['specific_internal'])
        native_change+=.5*mass*u[...,0]*np.sum(f['velocity']**2,axis=-1)
        physical_change=stable_internal_change(n,np.stack([sample(u[...,k]) for k in (4,5)],axis=-1),initial_n,initial_s,E*initial_t/(gamma-1))+kinetic
        momentum=float(mass*np.sum(fluid.volume*u[...,3]));expected_momentum=(1-(-1)**mode)*mass*density*np.pi*radius**2*amplitude*speed/wave*np.sin(omega*when)
        return dict(time=when,density_error_over_perturbation=norm(n-expected_n)/(amplitude*density),
            velocity_error_over_perturbation=norm(v-expected_v)/(amplitude*speed),
            pressure_error_over_perturbation=norm(p-expected_p)/(gamma*amplitude*pressure),
            particle_relative_change=float(np.sum(fluid.volume*(u[...,0]-initial_state[...,0]))/(density*volume)),
            native_energy_j=float(np.sum(fluid.volume*f['energy'])),physical_energy_j=float(volume*np.sum(weights*(kinetic+internal))),
            native_energy_change_stable_j=float(np.sum(fluid.volume*native_change)),physical_energy_change_stable_j=float(volume*np.sum(weights*physical_change)),
            axial_momentum_kg_m_s=momentum,linear_axial_momentum_kg_m_s=float(expected_momentum),
            momentum_error_over_perturbation=float(abs(momentum-expected_momentum)/(mass*density*volume*amplitude*speed)),
            end_pressure_force_n=diagnostic['end_pressure_force_n'],end_constraint_force_n=diagnostic['end_constraint_force_n'],
            axial_momentum_rate_n=diagnostic['axial_momentum_rate_n'],
            boundary_momentum_balance_n=diagnostic['axial_momentum_rate_n']-diagnostic['end_pressure_force_n']-diagnostic['end_constraint_force_n'],
            native_instantaneous_power_w=float(np.sum(fluid.volume*diagnostic['physical_energy_rate'])),
            minimum_density_m3=float(u[...,0].min()),maximum_axial_velocity_m_s=float(abs(f['velocity'][...,2]).max()),
            maximum_transverse_velocity_m_s=float(abs(f['velocity'][...,:2]).max()),
            axial_profile=dict(z_m=z.tolist(),density_m3=n.tolist(),velocity_m_s=v.tolist(),pressure_pa=p.tolist()))
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:
            raise RuntimeError('Finite-tube run exceeded own4GiB/256MiB resource limits')
    np.savez_compressed(destination/'initial.npz',state=state);record['initial_sha256']=sha(destination/'initial.npz')
    initial=observe(state,0.);record['initial']=initial;save();when=0.
    # Use the acoustic perturbation energy scale, not the large thermal background.
    energy_scale=.25*mass*density*volume*(amplitude*speed)**2
    record['acoustic_energy_scale_j']=energy_scale
    stride=max(1,steps//80)
    for step in range(1,steps+1):
        try:updated=fluid.step(state,dt);resources()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error));break
        when=step*dt;state=updated;raw=destination/f'step_{step:05}.npz';np.savez_compressed(raw,state=state)
        record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            measured=observe(state,when);measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),
                native_energy_change_over_acoustic_scale=measured['native_energy_change_stable_j']/energy_scale,
                physical_energy_change_over_acoustic_scale=measured['physical_energy_change_stable_j']/energy_scale)
            record['records'].append(measured);save()
        if step==1 or step%(steps//8 or 1)==0:print('finite tube',axial_points,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=when)
    if provenance(['examples.finite_cylinder_acoustic'])!=identity:raise RuntimeError('Source changed during finite-tube run')
    save();print('finite tube terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--axial-points',type=int,default=12)
    parser.add_argument('--mode',type=int,default=1);parser.add_argument('--steps',type=int,default=1600);parser.add_argument('--amplitude',type=float,default=1e-5)
    args=parser.parse_args();result=run(args.destination,axial_points=args.axial_points,steps=args.steps,amplitude=args.amplitude,mode=args.mode)
    raise SystemExit(0 if result['success'] else 1)
