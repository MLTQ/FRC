"""Complete heated material/current expansion with physical species transport."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,electron_mass as ME,mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_current_space import MaterialCurrentSpace
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_thermal_current import MaterialThermalCurrent
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,*,nr=4,nz=8,degree=3,quadrature=6,steps=2000,resistivity=.001,transport='upwind'):
    out=Path(destination)
    if out.exists() or degree not in (1,3,5) or min(nr,nz)<degree or quadrature<2*degree or steps<4:
        raise ValueError('Fresh output, supported complete mesh/quadrature/time resolution required')
    started=time.monotonic();identity=provenance(['examples.material_thermal_expansion']);out.mkdir(parents=True);freeze(out,identity)
    radius=.035;length=.07;mi=39.948*atomic_mass;mass=mi+ME;nc=1e18;amplitude=2e-4
    mesh=MHDMesh(radius,length,nr,nz);domain=MaterialCylinder(mesh,quadrature=quadrature);space=MaterialCurrentSpace(mesh,degree=degree,quadrature=quadrature)
    r=space.reference_r.reshape(nr*nz,-1);z=space.reference_z.reshape(nr*nz,-1)
    hr=1-(r/radius)**2;hz=1-(z/length)**2;profile=hr*hz
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,nc*profile**1.5,np.broadcast_to(profile[...,None],profile.shape+(2,)),mi)
    fluid.initial.entropy[:]=-(2/3)*np.log(nc)
    model=MaterialThermalCurrent(fluid,space,resistivity=resistivity,coulomb_exchange=True,transport=transport)
    exact=amplitude*r/radius*profile**3
    current=np.stack((6*amplitude*r*z/(MU0*radius*length**2)*hr**3*hz**2,
        amplitude/(MU0*radius)*hz**3*(2*hr**3-6*(r/radius)**2*hr**2)),axis=-1)
    state=model.project_initial(fluid.initial,exact,current);initial=model.response(state)
    duration=radius/np.sqrt(10*E/mass);dt=duration/steps;energy0=initial['total_energy']
    inertia0=float(np.sum(mass*fluid.particles*np.sum(initial['fields']['points']**2,axis=-1)))
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,mode='free',nr=nr,nz=nz,degree=degree,quadrature=quadrature,
        radius_m=radius,half_length_m=length,ion_mass_kg=mi,central_density_m3=nc,amplitude_t=amplitude,source_sha256=identity['source_sha256'],records=[],
        initial_total_energy_j=energy0,initial_projection='magnetic_current',thermal_current=True,resistivity_ohm_m=resistivity,
        coulomb_exchange=True,electron_heat_fraction=model.electron_heat_fraction,entropy_transport=transport,
        scope='Simultaneous free material/current/species entropy with Hall, physical inertia/stress, pressure work, prescribed scalar resistance and its heat, conditional Maxwellian exchange. Toroidal B, exact empty exterior; not a validated heated edge, Coulomb conductivity, coil model or FRC formation.')
    wb_integral=0.;twb_integral=0.;heat_integral=0.;exchange_integral=0.;previous=initial
    def observe(u,when,a):
        fields=a['fields'];n=fields['n'].ravel();j=E*n[:,None]*a['drift']
        inertia=float(np.sum(mass*fluid.particles*np.sum(fields['points']**2,axis=-1)))
        predicted=inertia0+2*energy0*when**2-2*(when*wb_integral-twb_integral)
        entropy=np.sum(fluid.particles[...,None]*u.fluid.entropy/(fluid.gammas-1),axis=(0,1))
        return dict(time=when,total_energy_j=a['total_energy'],relative_total_energy_change=a['total_energy']/energy0-1,
            magnetic_energy_j=a['magnetic_energy'],current_energy_j=a['current_energy'],bulk_kinetic_energy_j=fields['kinetic_energy'],internal_energy_j=fields['internal_energy'].tolist(),
            joule_power_w=a['joule_power'],joule_heat_j=heat_integral,exchange_power_to_electrons_w=a['exchange_power_to_electrons'],exchange_heat_to_electrons_j=exchange_integral,
            species_entropy_over_kB=entropy.tolist(),entropy_rates_over_kB=a['entropy_rates_over_kB'].tolist(),
            minimum_specific_entropy=np.min(u.fluid.entropy,axis=(0,1)).tolist(),maximum_transport_outflow_s=a['maximum_transport_outflow_s'],
            peak_b_t=float(abs(a['magnetic_field']).max()),peak_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),peak_drift_m_s=float(np.linalg.norm(a['drift'],axis=-1).max()),
            minimum_density_m3=float(n.min()),minimum_temperature_ev=float(fields['temperature'].min()),
            particle_mean_temperature_ev=(np.sum(fluid.particles[...,None]*fields['temperature'],axis=(0,1))/fluid.particles.sum()).tolist(),
            peak_temperature_ev=np.max(fields['temperature'],axis=(0,1)).tolist(),minimum_jacobian_m2=float(fields['jacobian'].min()),
            virial_inertia_error_over_initial=(inertia-predicted)/inertia0,residuals=a['residuals'],weak_charge_residual=a['weak_charge_residual'],pair_power_residual_w=a['pair_power_residual'],
            exchange_frequency_range_s=a['exchange_frequency_range_s'],exchange_coulomb_log_range=a['exchange_coulomb_log_range'],boundary_coordinates=u.fluid.coordinates[domain.boundary].tolist())
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def arrays(u):return dict(magnetic=u.magnetic,**{k:getattr(u.fluid,k) for k in u.fluid.__dataclass_fields__})
    np.savez_compressed(out/'initial.npz',**arrays(state),particles=fluid.particles,ids=fluid.ids)
    record['initial']=observe(state,0.,initial);save();stride=max(1,steps//80)
    for step in range(1,steps+1):
        t0=(step-1)*dt
        try:
            updated=model.step(state,dt);response=model.response(updated)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:
                raise RuntimeError('Own4GiB/256MiB resource guard')
        except Exception as exc:
            record.update(failed_step=step,failure_type=type(exc).__name__,failure_message=str(exc));break
        state=updated;when=step*dt
        wb_integral+=.5*dt*(previous['magnetic_energy']+response['magnetic_energy'])
        twb_integral+=.5*dt*(t0*previous['magnetic_energy']+when*response['magnetic_energy'])
        heat_integral+=.5*dt*(previous['joule_power']+response['joule_power'])
        exchange_integral+=.5*dt*(previous['exchange_power_to_electrons']+response['exchange_power_to_electrons']);previous=response
        raw=out/f'step_{step:05}.npz';np.savez_compressed(raw,**arrays(state));record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            row=observe(state,when,response);row.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(row);save()
        if step==1 or step%(steps//8 or 1)==0:print('material thermal',nr,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.material_thermal_expansion'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material thermal terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--nz',type=int,default=8)
    p.add_argument('--degree',type=int,default=3);p.add_argument('--quadrature',type=int,default=6);p.add_argument('--steps',type=int,default=2000);p.add_argument('--resistivity',type=float,default=.001);p.add_argument('--transport',choices=('skew','upwind'),default='upwind')
    a=p.parse_args();d=run(a.destination,nr=a.nr,nz=a.nz,degree=a.degree,quadrature=a.quadrature,steps=a.steps,resistivity=a.resistivity,transport=a.transport)
    raise SystemExit(0 if d['success'] else 1)
