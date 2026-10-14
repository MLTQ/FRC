"""Complete ideal magnetic free expansion; independent initial fields and virial."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,electron_mass as ME,mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_magnetic import CompactMagneticField,MaterialMagnetic,MaterialIdealPlasma
from .implicit_plasma_records import provenance,freeze,sha


def spherical_fields(r,z,radius,bp,bt):
    h=1-(r*r+z*z)/radius**2
    b=np.stack((4*bp*r*z/radius**2*h**3,bt*r/radius*h**3,bp*h**4-4*bp*r*r/radius**2*h**3),axis=-1)
    j=np.stack((6*bt*r*z/radius**3*h**2,4*bp*r/radius**2*h**2*(5-11*(r*r+z*z)/radius**2),
        bt/radius*(2*h**3-6*r*r/radius**2*h**2)),axis=-1)/MU0
    return b,j


def run(destination,*,nr=8,steps=800):
    destination=Path(destination)
    if destination.exists() or nr<2 or steps<4:raise ValueError('Fresh output and resolved positive trajectory required')
    started=time.monotonic();identity=provenance(['examples.material_magnetic_expansion']);destination.mkdir(parents=True);freeze(destination,identity)
    radius=.035;density=1e24;temperature=np.array([1.,1.]);mi=39.948*atomic_mass;mass=mi+ME;bp=.3;bt=.15
    coefficient=5*E*temperature.sum()/mass;speed=np.sqrt(coefficient);duration=radius/speed;dt=duration/steps
    domain=MovingDomain(nr,2*nr,radius=radius,matching_radius=4*radius,vacuum_layers=2,modes=1,quadrature=6)
    q,_=np.polynomial.legendre.leggauss(domain.quadrature);xi,tau=np.meshgrid(.5*(q+1),.5*(q+1),indexing='ij');s=(np.repeat(np.arange(nr),domain.nt)[:,None]+xi.ravel()[None,:])/nr
    h=1-s*s;fluid=MaterialTwoSpecies(domain,domain.fluid_reference,density*h**1.5,h[...,None]*temperature,mi)
    magnetic=MaterialMagnetic(fluid,CompactMagneticField(fluid,poloidal_tesla=bp,toroidal_tesla=bt));model=MaterialIdealPlasma(fluid,magnetic)
    state=model.initial;initial=model.fields(state);energy0=initial['total_energy'];particles=float(fluid.particles.sum())
    def inertia(f):return float(np.sum(mass*fluid.particles*np.sum(f['points']**2,axis=-1)))
    moment0=inertia(initial);moment_rate0=2*np.sum(mass*fluid.particles*np.sum(initial['points']*initial['sampled_velocity'],axis=-1))
    # Independent continuum spherical quadrature, unrelated to production cells.
    aq,aw=np.polynomial.legendre.leggauss(80);rad=.5*(aq+1)*radius;rw=.5*aw*radius
    rr=rad[:,None]*np.sqrt(1-aq[None,:]**2);zz=rad[:,None]*aq[None,:]
    exact_b,exact_j=spherical_fields(rr,zz,radius,bp,bt);weight=2*np.pi*rad[:,None]**2*rw[:,None]*aw[None,:]
    exact_magnetic=float(np.sum(weight*np.sum(exact_b**2,axis=-1))/(2*MU0))
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,radial_cells=nr,angular_cells=domain.nt,
        quadrature=domain.quadrature,radius_m=radius,central_density_m3=density,central_temperature_ev=temperature.tolist(),ion_mass_kg=mi,
        poloidal_amplitude_t=bp,toroidal_amplitude_t=bt,speed_scale_m_s=speed,source_sha256=identity['source_sha256'],records=[],
        native_initial_energy_j=energy0,exact_initial_magnetic_energy_j=exact_magnetic,native_initial_magnetic_energy_j=initial['magnetic_energy'],
        relative_initial_magnetic_energy_error=initial['magnetic_energy']/exact_magnetic-1,
        continuum='Formal ideal MHD: full magnetic deformation/torque, adiabatic species carried with bulk; no finite-current inertia, Hall, relative species transport, resistivity or external coils. Mathematical compact initial field, no exterior matter/field. No formation or laboratory validity claim.')
    integral_b=0.;integral_tb=0.;previous_b=initial['magnetic_energy']
    def observe(u,t,f):
        physical=magnetic.samples(u);b=physical['magnetic_field'];j=physical['current_density'];vol=f['volume'];volume=vol.sum()
        w=j/(E*f['n'][...,None]);drift=np.linalg.norm(w,axis=-1)
        current_energy=.5*mi*ME/mass*np.sum(fluid.particles*np.sum(w*w,axis=-1))
        ion_skin=np.sqrt(mi/(MU0*f['n']*E*E));moment=inertia(f)
        moment_reference=moment0+moment_rate0*t+2*energy0*t*t-2*(t*integral_b-integral_tb)
        virial_rate=4*f['kinetic_energy']+4*f['internal_energy'].sum()+2*f['magnetic_energy']
        return dict(time=t,total_energy_j=f['total_energy'],relative_energy_change=f['total_energy']/energy0-1,
            magnetic_energy_j=f['magnetic_energy'],kinetic_energy_j=f['kinetic_energy'],internal_energy_j=f['internal_energy'].tolist(),
            inertia_kg_m2=moment,virial_inertia_reference_kg_m2=float(moment_reference),
            virial_inertia_error_over_initial=(moment-moment_reference)/moment0,virial_second_rate_j=float(virial_rate),
            peak_b_t=float(np.linalg.norm(b,axis=-1).max()),peak_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),
            peak_speed_m_s=float(np.sqrt(np.sum(f['sampled_velocity']**2,axis=-1)+(f['points'][...,0]*f['sampled_omega'])**2).max()),
            peak_toroidal_speed_m_s=float(abs(f['points'][...,0]*f['sampled_omega']).max()),
            maximum_div_b_over_b_radius=float(abs(physical['magnetic_divergence']).max()*radius/bp),
            maximum_reference_pullback_difference_t=float(abs(b-f['magnetic_field']).max()),
            angular_momentum_kg_m2_s=float(u.angular_momentum.sum()),axial_momentum_kg_m_s=float(u.momentum[:,1].sum()),
            minimum_jacobian_m2=float(f['jacobian'].min()),minimum_density_m3=float(f['n'].min()),minimum_temperature_ev=float(f['temperature'].min()),
            inferred_current_kinetic_energy_j=float(current_energy),inferred_current_energy_over_total=float(current_energy/energy0),
            ion_skin_depth_range_over_initial_radius=[float(ion_skin.min()/radius),float(ion_skin.max()/radius)],
            maximum_drift_over_ion_thermal_speed=float(np.max(drift/np.sqrt(E*f['temperature'][...,1]/mi))),
            magnetic_rms_t=float(np.sqrt(np.sum(vol*np.sum(b*b,axis=-1))/volume)),
            current_rms_a_m2=float(np.sqrt(np.sum(vol*np.sum(j*j,axis=-1))/volume)),
            cell_current_rms_a_m2=float(np.sqrt(np.sum(vol*np.sum(physical['cell_current_density']**2,axis=-1))/volume)),
            projected_minus_cell_current_rms_a_m2=float(np.sqrt(np.sum(vol*np.sum((j-physical['cell_current_density'])**2,axis=-1))/volume)),
            boundary_coordinates=u.coordinates[domain.boundary].tolist())
    x=initial['points'];analytic_b,analytic_j=spherical_fields(x[...,0],x[...,1],radius,bp,bt);sampled=magnetic.samples(state);v0=initial['volume']
    record['initial_b_rms_error_over_amplitude']=float(np.sqrt(np.sum(v0*np.sum((sampled['magnetic_field']-analytic_b)**2,axis=-1))/v0.sum())/bp)
    record['initial_cell_current_rms_error_over_scale']=float(np.sqrt(np.sum(v0*np.sum((sampled['cell_current_density']-analytic_j)**2,axis=-1))/v0.sum())/(bp/(MU0*radius)))
    record['initial_current_rms_error_over_scale']=float(np.sqrt(np.sum(v0*np.sum((sampled['current_density']-analytic_j)**2,axis=-1))/v0.sum())/(bp/(MU0*radius)))
    # Pointwise continuum acceleration includes independent pressure and JxB.
    exact_h=1-np.sum(x*x,axis=-1)/radius**2;exact_density=density*exact_h**1.5
    acceleration=np.cross(analytic_j,analytic_b)/(mass*exact_density[...,None]);acceleration[...,[0,2]]+=coefficient*x/radius**2
    rhs=model.derivative(state);nodal=np.zeros_like(state.coordinates)
    nodal[fluid.active_r,0]=fluid.mass_factors[0].solve(rhs.momentum[fluid.active_r,0]);nodal[:,1]=fluid.mass_factors[1].solve(rhs.momentum[:,1])
    sampled_acceleration=np.zeros_like(acceleration);sampled_acceleration[...,[0,2]]=np.einsum('qn,cna->cqa',fluid.shape,nodal[fluid.ids])
    from scipy.sparse.linalg import splu
    angular_acceleration=splu(fluid.gram(mass*fluid.particles*x[...,0]**2)).solve(rhs.angular_momentum)
    sampled_acceleration[...,1]=x[...,0]*np.einsum('qn,cn->cq',fluid.shape,angular_acceleration[fluid.ids])
    record['initial_acceleration_particle_rms_error_over_scale']=float(np.sqrt(np.sum(fluid.particles*np.sum((sampled_acceleration-acceleration)**2,axis=-1))/particles)/(coefficient/radius))
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    np.savez_compressed(destination/'initial.npz',**{k:getattr(state,k) for k in state.__dataclass_fields__},particles=fluid.particles,ids=fluid.ids,profile_s=s)
    record['initial_sha256']=sha(destination/'initial.npz');record['initial']=observe(state,0.,initial);save();stride=max(1,steps//80)
    for step in range(1,steps+1):
        try:
            updated=model.step(state,dt);f=model.fields(updated)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
        except Exception as exc:
            record.update(failed_step=step,failure_type=type(exc).__name__,failure_message=str(exc));break
        state=updated;when=step*dt;now_b=f['magnetic_energy']
        integral_b+=.5*dt*(previous_b+now_b);integral_tb+=.5*dt*((when-dt)*previous_b+when*now_b);previous_b=now_b
        raw=destination/f'step_{step:05}.npz';np.savez_compressed(raw,**{k:getattr(state,k) for k in state.__dataclass_fields__});record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            measured=observe(state,when,f);measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(measured);save()
        if step==1 or step%(steps//8 or 1)==0:print('material magnetic',nr,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.material_magnetic_expansion'])!=identity:raise RuntimeError('Source changed during trajectory')
    save();print('material magnetic terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=8);p.add_argument('--steps',type=int,default=800);a=p.parse_args()
    d=run(a.destination,nr=a.nr,steps=a.steps);raise SystemExit(0 if d['success'] else 1)
