"""Complete free-surface two-species trajectories against affine Euler solutions."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,electron_mass as ME
from scipy.integrate import solve_ivp
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.material_two_species import MaterialTwoSpecies
from .implicit_plasma_records import provenance,freeze,sha


def reference(when,*,radius,length,coefficient,rotation):
    """Independent affine Euler axes; sphere/no-spin has an exact closed form."""
    if length==radius and rotation==0:
        a=np.sqrt(radius**2+coefficient*when**2);ad=coefficient*when/a
        return np.array([a,a,ad,ad,0.,0.])
    def rhs(t,y):
        a,c,ad,cd,phase=y;volume=a*a*c/(radius*radius*length);thermal=coefficient*volume**(-2/3)
        return [ad,cd,thermal/a+rotation**2*radius**4/a**3,thermal/c,rotation*(radius/a)**2]
    if when==0:return np.array([radius,length,0.,0.,rotation,0.])
    result=solve_ivp(rhs,(0,when),[radius,length,0.,0.,0.],method='DOP853',rtol=2e-13,atol=1e-14)
    if not result.success:raise RuntimeError(result.message)
    a,c,ad,cd,phase=result.y[:,-1];return np.array([a,c,ad,cd,rotation*(radius/a)**2,phase])


def run(destination,*,nr=8,nt=None,steps=800,mode='sphere'):
    destination=Path(destination);nt=2*nr if nt is None else nt
    if destination.exists() or nr<2 or nt<4 or steps<4 or mode not in ('sphere','ellipsoid','rotating'):
        raise ValueError('Fresh output, resolved material mesh and complete positive trajectory required')
    started=time.monotonic();identity=provenance(['examples.material_vacuum_expansion']);destination.mkdir(parents=True);freeze(destination,identity)
    radius=.035;length=radius if mode=='sphere' else .07;density=1e20;temperature=np.array([1.,1.]);mi=39.948*atomic_mass;mass=mi+ME
    coefficient=5*E*temperature.sum()/mass;duration=radius/np.sqrt(coefficient);dt=duration/steps
    rotation=.5*np.sqrt(coefficient)/radius if mode=='rotating' else 0.
    domain=MovingDomain(nr,nt,radius=radius,matching_radius=4*length,vacuum_layers=2,modes=1,quadrature=6)
    coordinates=domain.fluid_reference.copy();coordinates[:,1]*=length/radius
    q,_=np.polynomial.legendre.leggauss(domain.quadrature);xi,tau=np.meshgrid(.5*(q+1),.5*(q+1),indexing='ij');xi=xi.ravel()
    label_s=np.repeat(np.arange(nr),nt)[:,None]/nr+xi[None,:]/nr
    profile=1-label_s**2;n=density*profile**1.5;temperature_q=profile[...,None]*temperature
    fluid=MaterialTwoSpecies(domain,coordinates,n,temperature_q,mi,angular_velocity=np.full(domain.fluid_node_count,rotation))
    state=fluid.initial;initial=fluid.fields(state);particles=float(np.sum(fluid.particles));exact_particles=np.pi**2/8*radius**2*length*density
    energy0=initial['energy'];exact_energy=exact_particles*E*temperature.sum()*15/16+mass*exact_particles*rotation**2*radius**2/8
    initial_angular=float(state.angular_momentum.sum());speed_scale=np.sqrt(coefficient)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,mode=mode,
        radial_cells=nr,angular_cells=nt,quadrature=domain.quadrature,initial_radial_axis_m=radius,initial_axial_axis_m=length,
        central_density_m3=density,central_temperature_ev=temperature.tolist(),ion_mass_kg=mi,rotation_s=rotation,
        speed_scale_m_s=speed_scale,source_sha256=identity['source_sha256'],records=[],
        exact_particles=exact_particles,particle_quadrature=particles,relative_initial_particle_error=(particles-exact_particles)/exact_particles,
        exact_initial_energy_j=exact_energy,native_initial_energy_j=energy0,relative_initial_energy_error=(energy0-exact_energy)/exact_energy,
        continuum='gamma_e=gamma_i=5/3; physical total particle mass mi+me; zero B/current, equal electron/ion temperature, so Hall/Joule/current-stress/exchange terms vanish in this mathematical fixture. No finite-field transport or formation claim.',
        reference='Independent physical-vacuum affine Euler: n0=nc*(1-s^2)^(3/2),T0=Tc*(1-s^2); a_ddot=C*D^(-2/3)/a+Omega0^2*a0^4/a^3,c_ddot=C*D^(-2/3)/c; C=5e*(Te0+Ti0)/(mi+me),D=a^2*c/(a0^2*c0). Sphere: a=c=sqrt(R0^2+C*t^2).')
    def observe(u,t):
        f=fluid.fields(u);a,c,ad,cd,omega,phase=reference(t,radius=radius,length=length,coefficient=coefficient,rotation=rotation)
        ratio=a*a*c/(radius*radius*length)
        expected_nodes=coordinates*np.array([a/radius,c/length]);expected_velocity=u.coordinates*np.array([ad/a,cd/c])
        expected_q_velocity=f['points']*np.array([ad/a,cd/c]);q_velocity_error=np.sum((f['sampled_velocity']-expected_q_velocity)**2,axis=-1)
        expected_density=n/ratio;expected_temperature=temperature_q*ratio**(-2/3)
        # This evaluates the independent physical-space compact support. The
        # positive part belongs only to the analytic reference, never state repair.
        h=np.maximum(1-(f['points'][...,0]/a)**2-(f['points'][...,1]/c)**2,0.)
        euler_density=density/ratio*h**1.5;euler_temperature=h[...,None]*temperature*ratio**(-2/3)
        weight=f['volume'];volume=float(np.sum(weight));norm=lambda v:float(np.sqrt(np.sum(weight*v*v)/volume))
        boundary=u.coordinates[domain.boundary];boundary_reference=expected_nodes[domain.boundary]
        qvr=f['sampled_velocity'];qomega=f['sampled_omega']
        expected_kinetic=mass*exact_particles/16*(2*ad**2+cd**2+2*a*a*omega*omega)
        expected_internal=exact_particles*E*temperature.sum()*15/16*ratio**(-2/3)
        return dict(time=t,radial_axis_reference_m=a,axial_axis_reference_m=c,
            maximum_boundary_position_error_over_initial_radius=float(np.max(np.linalg.norm(boundary-boundary_reference,axis=1))/radius),
            material_position_rms_error_over_initial_radius=float(np.sqrt(np.sum(fluid.particles*np.sum((f['points']-np.einsum('qn,cna->cqa',fluid.shape,expected_nodes[fluid.ids]))**2,axis=-1))/particles)/radius),
            velocity_rms_error_over_speed=float(np.sqrt(np.sum(fluid.particles*q_velocity_error)/particles)/speed_scale),
            material_angle_rms_error_rad=float(np.sqrt(np.sum(fluid.particles*(np.einsum('qn,cn->cq',fluid.shape,u.angle[fluid.ids])-phase)**2)/particles)),
            reference_angle_rad=float(phase),
            angular_velocity_rms_error_over_scale=float(np.sqrt(np.sum(fluid.particles*(qomega-omega)**2)/particles)/(speed_scale/radius)),
            material_density_rms_error_over_center=norm(f['n']-expected_density)/(density/ratio),
            spatial_density_rms_error_over_center=norm(f['n']-euler_density)/(density/ratio),
            material_temperature_rms_error_over_center=[norm(f['temperature'][...,k]-expected_temperature[...,k])/(temperature[k]*ratio**(-2/3)) for k in range(2)],
            spatial_temperature_rms_error_over_center=[norm(f['temperature'][...,k]-euler_temperature[...,k])/(temperature[k]*ratio**(-2/3)) for k in range(2)],
            relative_energy_change=(f['energy']-energy0)/energy0,kinetic_energy_j=f['kinetic_energy'],internal_energy_j=f['internal_energy'].tolist(),
            exact_kinetic_energy_j=float(expected_kinetic),exact_internal_energy_j=float(expected_internal),
            axial_momentum_kg_m_s=float(u.momentum[:,1].sum()),angular_momentum_kg_m2_s=float(u.angular_momentum.sum()),
            angular_momentum_change=float(u.angular_momentum.sum()-initial_angular),
            particle_reconstructed=float(np.sum(weight*f['n'])),minimum_interior_density_m3=float(f['n'].min()),
            minimum_interior_temperature_ev=float(f['temperature'].min()),minimum_jacobian_m2=float(f['jacobian'].min()),
            maximum_speed_m_s=float(np.sqrt(np.sum(qvr*qvr,axis=-1)+(f['points'][...,0]*qomega)**2).max()),
            boundary_coordinates=boundary.tolist())
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def guard():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Material expansion exceeded own4GiB/256MiB resource bounds')
    np.savez_compressed(destination/'initial.npz',**{k:getattr(state,k) for k in state.__dataclass_fields__},particles=fluid.particles,ids=fluid.ids,profile_s=label_s)
    record['initial_sha256']=sha(destination/'initial.npz');record['initial']=observe(state,0.);save();stride=max(1,steps//80)
    for step in range(1,steps+1):
        try:updated=fluid.step(state,dt);guard()
        except Exception as exc:
            record.update(failed_step=step,failure_type=type(exc).__name__,failure_message=str(exc));break
        state=updated;when=step*dt;raw=destination/f'step_{step:05}.npz';np.savez_compressed(raw,**{k:getattr(state,k) for k in state.__dataclass_fields__})
        record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            measured=observe(state,when);measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(measured);save()
        if step==1 or step%(steps//8 or 1)==0:print('material expansion',mode,nr,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.material_vacuum_expansion'])!=identity:raise RuntimeError('Material expansion source changed during run')
    save();print('material expansion terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=8);p.add_argument('--nt',type=int);p.add_argument('--steps',type=int,default=800);p.add_argument('--mode',choices=('sphere','ellipsoid','rotating'),default='sphere');a=p.parse_args()
    d=run(a.destination,nr=a.nr,nt=a.nt,steps=a.steps,mode=a.mode);raise SystemExit(0 if d['success'] else 1)
