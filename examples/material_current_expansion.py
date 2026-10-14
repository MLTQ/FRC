"""Independent inertial expansion and actual free material/current trajectories."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,electron_mass as ME,mu_0 as MU0
from scipy.special import jv,jn_zeros
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_current_space import MaterialCurrentSpace
from plasma_bench.material_two_species import MaterialTwoSpecies,MaterialSpeciesState
from plasma_bench.material_carrier_current import MaterialCarrierCurrent,MaterialCurrentState,positive_solve
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,*,nr=4,nz=8,steps=400,mode='reference',initial_projection='magnetic_current',degree=3,quadrature=6):
    destination=Path(destination)
    if destination.exists() or mode not in ('reference','free') or initial_projection not in ('magnetic','magnetic_current') or steps<4:raise ValueError('Fresh output and complete declared mode required')
    if degree not in (1,3,5) or min(nr,nz)<degree or quadrature<2*degree:raise ValueError('Supported degree, enough cells and quadrature >= 2*degree required')
    started=time.monotonic();identity=provenance(['examples.material_current_expansion']);destination.mkdir(parents=True);freeze(destination,identity)
    radius=.035;length=.07;mi=39.948*atomic_mass;mass=mi+ME;mu=mi*ME/mass
    mesh=MHDMesh(radius,length,nr,nz);domain=MaterialCylinder(mesh,quadrature=quadrature);space=MaterialCurrentSpace(mesh,degree=degree,quadrature=quadrature)
    r=space.reference_r.reshape(nr*nz,-1);z=space.reference_z.reshape(nr*nz,-1)
    kr=jn_zeros(1,1)[0]/radius;kz=np.pi/(2*length);density=mu*(kr*kr+kz*kz)/(MU0*E*E) if mode=='reference' else 1e18
    h=np.ones_like(r) if mode=='reference' else (1-(r/radius)**2)*(1-(z/length)**2)
    n=density*h**1.5;temperature=np.broadcast_to(h[...,None],h.shape+(2,)).copy()
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,n,temperature,mi)
    # These supplied data are isentropic by definition; represent their exact
    # common entropy directly rather than retain exp/log construction roundoff.
    fluid.initial.entropy[:]=-(2/3)*np.log(density)
    model=MaterialCarrierCurrent(fluid,space);duration=1e-6 if mode=='reference' else radius/np.sqrt(10*E/mass);dt=duration/steps
    amplitude=1e-12 if mode=='reference' else 2e-4
    exact0=amplitude*jv(1,kr*r)*np.cos(kz*z) if mode=='reference' else amplitude*r/radius*h**3
    if mode=='reference':
        reference_jr=amplitude*kz/MU0*jv(1,kr*r)*np.sin(kz*z)
        reference_jz=amplitude*kr/MU0*jv(0,kr*r)*np.cos(kz*z)
    else:
        hr=1-(r/radius)**2;hz=1-(z/length)**2
        reference_jr=6*amplitude*r*z/(MU0*radius*length**2)*hr**3*hz**2
        reference_jz=amplitude/(MU0*radius)*hz**3*(2*hr**3-6*(r/radius)**2*hr**2)
    current0=np.stack((reference_jr,reference_jz),axis=-1)
    if initial_projection=='magnetic_current':state=model.project_initial(fluid.initial,exact0,current0)
    else:
        geom=model.geometry(fluid.initial);b=positive_solve(geom['h'],np.asarray(geom['bmap'].T@(geom['volume']/MU0*exact0.ravel())))
        state=MaterialCurrentState(fluid.initial,b)
    b=state.magnetic;x0=fluid.initial.coordinates.copy()
    def prescribed(when,coeff):
        lam=1+when/duration;velocity=x0/duration
        momentum=np.column_stack((fluid.mass_matrix@velocity[:,0],fluid.mass_matrix@velocity[:,1]));momentum[domain.axis,0]=0
        return MaterialCurrentState(MaterialSpeciesState(lam*x0,momentum,np.zeros(fluid.nodes),np.zeros(fluid.nodes),fluid.initial.entropy.copy()),coeff)
    if mode=='reference':state=prescribed(0.,b)
    initial=model.response(state);energy0=initial['total_energy'];field_energy0=initial['magnetic_energy']+initial['current_energy']
    inertia0=float(np.sum(mass*fluid.particles*np.sum(initial['fields']['points']**2,axis=-1)))
    inertia_rate0=float(2*np.sum(mass*fluid.particles*np.sum(initial['fields']['points']*initial['fields']['sampled_velocity'],axis=-1)))
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,mode=mode,nr=nr,nz=nz,degree=degree,quadrature=quadrature,
        radius_m=radius,half_length_m=length,ion_mass_kg=mi,central_density_m3=density,amplitude_t=amplitude,source_sha256=identity['source_sha256'],records=[],
        initial_field_energy_j=field_energy0,initial_total_energy_j=energy0,initial_projection=initial_projection,
        scope='Moving finite-mass Hall/current-inertia/stress sector; Btheta and poloidal fluid/carrier velocity. Equal constant species entropies, eta=0. Reference motion prescribed; free mode evolves fluid simultaneously. No exterior matter, coil or FRC formation claim.')
    init_volume=initial['volume'];init_current=E*initial['fields']['n'].ravel()[:,None]*initial['drift'];target_current=current0.reshape(-1,2)
    record['initial_magnetic_relative_l2_error']=float(np.sqrt(np.sum(init_volume*(initial['magnetic_field']-exact0.ravel())**2)/np.sum(init_volume*exact0.ravel()**2)))
    record['initial_current_relative_l2_error']=float(np.sqrt(np.sum(init_volume[:,None]*(init_current-target_current)**2)/np.sum(init_volume[:,None]*target_current**2)))
    work=0.;wb_integral=0.;twb_integral=0.;prev=initial
    def observe(u,t,a):
        n=a['fields']['n'].ravel();j=E*n[:,None]*a['drift'];v=a['fields']['volume'].ravel();volume=v.sum()
        lam=1+t/duration;decay=2/(1+lam)
        expected_b=exact0.ravel()*decay/lam**2;expected_j=np.column_stack((reference_jr.ravel(),reference_jz.ravel()))*decay/lam**3
        rms=lambda val:float(np.sqrt(np.sum(v*val)/volume))
        moment=float(np.sum(mass*fluid.particles*np.sum(a['fields']['points']**2,axis=-1)))
        imodel=inertia0+inertia_rate0*t+2*energy0*t*t-2*(t*wb_integral-twb_integral)
        row=dict(time=t,total_energy_j=a['total_energy'],relative_total_energy_change=a['total_energy']/energy0-1,
            magnetic_energy_j=a['magnetic_energy'],current_energy_j=a['current_energy'],bulk_kinetic_energy_j=a['fields']['kinetic_energy'],internal_energy_j=a['fields']['internal_energy'].tolist(),
            field_current_work_residual_over_initial=(a['magnetic_energy']+a['current_energy']-field_energy0+work)/field_energy0,
            pair_power_residual_w=a['pair_power_residual'],residuals=a['residuals'],weak_charge_residual=a['weak_charge_residual'],
            peak_b_t=float(abs(a['magnetic_field']).max()),peak_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),
            peak_drift_m_s=float(np.linalg.norm(a['drift'],axis=-1).max()),minimum_density_m3=float(n.min()),
            minimum_temperature_ev=float(a['fields']['temperature'].min()),minimum_jacobian_m2=float(a['fields']['jacobian'].min()),
            virial_inertia_error_over_initial=(moment-imodel)/inertia0,
            boundary_coordinates=u.fluid.coordinates[domain.boundary].tolist())
        if mode=='reference':row.update(magnetic_rms_error_over_initial_amplitude=rms((a['magnetic_field']-expected_b)**2)/amplitude,
            current_rms_error_over_initial_scale=rms(np.sum((j-expected_j)**2,axis=-1))/(amplitude*np.hypot(kr,kz)/MU0),
            expected_material_amplitude_ratio=decay,expected_physical_amplitude_ratio=decay/lam**2)
        return row
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def arrays(u):return dict(magnetic=u.magnetic,**{k:getattr(u.fluid,k) for k in u.fluid.__dataclass_fields__})
    np.savez_compressed(destination/'initial.npz',**arrays(state),particles=fluid.particles,ids=fluid.ids)
    record['initial']=observe(state,0.,initial);save();stride=max(1,steps//80)
    for step in range(1,steps+1):
        t0=(step-1)*dt
        try:
            if mode=='reference':
                rhs=lambda at,bb:model.response(prescribed(at,bb))['magnetic_rate']
                k1=rhs(t0,b);k2=rhs(t0+.5*dt,b+.5*dt*k1);k3=rhs(t0+.5*dt,b+.5*dt*k2);k4=rhs(t0+dt,b+dt*k3)
                updated=prescribed(step*dt,b+dt/6*(k1+2*k2+2*k3+k4))
            else:updated=model.step(state,dt)
            response=model.response(updated)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
        except Exception as exc:
            record.update(failed_step=step,failure_type=type(exc).__name__,failure_message=str(exc));break
        state=updated;b=state.magnetic;when=step*dt;work+=.5*dt*(prev['mechanical_power']+response['mechanical_power'])
        wb_integral+=.5*dt*(prev['magnetic_energy']+response['magnetic_energy']);twb_integral+=.5*dt*(t0*prev['magnetic_energy']+when*response['magnetic_energy']);prev=response
        raw=destination/f'step_{step:05}.npz';np.savez_compressed(raw,**arrays(state));record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            row=observe(state,when,response);row.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(row);save()
        if step==1 or step%(steps//8 or 1)==0:print('material current',mode,nr,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.material_current_expansion'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material current terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--nz',type=int,default=8);p.add_argument('--steps',type=int,default=400);p.add_argument('--mode',choices=('reference','free'),default='reference');p.add_argument('--initial-projection',choices=('magnetic','magnetic_current'),default='magnetic_current');p.add_argument('--degree',type=int,default=3);p.add_argument('--quadrature',type=int,default=6);a=p.parse_args()
    d=run(a.destination,nr=a.nr,nz=a.nz,steps=a.steps,mode=a.mode,initial_projection=a.initial_projection,degree=a.degree,quadrature=a.quadrature);raise SystemExit(0 if d['success'] else 1)
