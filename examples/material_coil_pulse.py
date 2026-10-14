"""A complete prescribed coil pulse through the actual freely moving plasma."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_magnetic_vacuum import MaterialMagneticVacuum
from plasma_bench.material_vacuum_space import MaterialVacuumSpace
from plasma_bench.material_coil_current import MaterialCoilCurrent
from plasma_bench.material_split_current import SplitMaterialCoilCurrent
from plasma_bench.material_carrier_current import MaterialCurrentState
from plasma_bench.electromagnetics import loop_field
from plasma_bench.waveforms import Unipolar,to_dict
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,*,nr=4,steps=2000,peak_field=1e-5,integrator='rk4'):
    out=Path(destination)
    if out.exists() or steps<4 or not np.isfinite(peak_field) or integrator not in ('rk4','split'):raise ValueError('Fresh complete pulse directory, supported integrator and finite amplitude required')
    started=time.monotonic();identity=provenance(['examples.material_coil_pulse']);out.mkdir(parents=True);freeze(out,identity)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,nr=nr,nz=2*nr,degree=3,quadrature=6,vacuum_cells=4,modes=6,
        source_sha256=identity['source_sha256'],radius_m=.035,half_length_m=.07,central_density_m3=1e18,resistivity_ohm_m=.001,
        nominal_vacuum_center_peak_t=peak_field,entropy_transport='upwind',integrator=integrator,records=[],
        scope='Prescribed finite circular coil pulse, all finite-mass moving plasma/current/species terms and exact matter-free magnetic exterior. Coarse exploratory model; no converged formation, physical transport/ionization, edge or laboratory-design claim.')
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        R=record['radius_m'];L=record['half_length_m'];nc=record['central_density_m3'];mesh=MHDMesh(R,L,nr,2*nr)
        domain=MaterialCylinder(mesh,quadrature=6);space=MaterialVacuumSpace(mesh,degree=3,quadrature=6)
        vacuum=MaterialMagneticVacuum(mesh,degree=3,quadrature_order=6,vacuum_cells=4,modes=6)
        r,z=space.reference_r,space.reference_z;profile=((1-(r/R)**2)*(1-(z/L)**2)).reshape(nr*2*nr,-1)
        fluid=MaterialTwoSpecies(domain,domain.fluid_reference,nc*profile**1.5,np.repeat(profile[...,None],2,axis=-1),39.948*atomic_mass)
        fluid.initial.entropy[:]=-(2/3)*np.log(nc)
        coils=((.07,-.035,1.),(.07,.035,1.));center_per_amp=sum(float(loop_field(0.,0.,radius,center)[2])*turns for radius,center,turns in coils)
        waveform=Unipolar(peak_field/center_per_amp,rise=1.5e-6,hold=1e-6,fall=2e-6)
        solver=SplitMaterialCoilCurrent if integrator=='split' else MaterialCoilCurrent
        model=solver(fluid,space,vacuum,coils,waveform,resistivity=.001)
        state=MaterialCurrentState(fluid.initial,np.zeros(space.magnetic_count));initial=model.response(state,time=0.)
        duration=R/np.sqrt(10*E/fluid.mass);dt=duration/steps;energy0=initial['total_energy'];external_work=0.;source_work=0.;heat=0.;previous=initial
        record.update(duration=duration,dt=dt,ion_mass_kg=fluid.ion_mass,coils=coils,waveform=to_dict(waveform),initial_total_energy_j=energy0)
        def observe(u,t,response):
            f=response['fields'];n=f['n'].ravel();j=E*n[:,None]*response['drift']
            return dict(time=t,source_current_a=response['source_current_a'],source_linkage_wb=response['source_linkage_wb'],mutual_energy_j=response['mutual_energy'],
                total_energy_j=response['total_energy'],response_energy_j=response['response_energy'],external_work_j=external_work,source_response_work_j=source_work if integrator=='rk4' else None,
                relative_external_work_residual=(response['total_energy']-energy0-external_work)/energy0,
                relative_response_work_residual=(response['response_energy']-energy0-source_work)/energy0 if integrator=='rk4' else None,
                magnetic_induced_energy_j=response['magnetic_energy'],vacuum_induced_energy_j=response['vacuum_magnetic_energy'],current_energy_j=response['current_energy'],
                bulk_kinetic_energy_j=f['kinetic_energy'],internal_energy_j=f['internal_energy'].tolist(),joule_heat_j=heat,
                peak_b_components_t=abs(response['magnetic_field']).max(axis=0).tolist(),peak_induced_b_components_t=abs(response['induced_magnetic_field']).max(axis=0).tolist(),
                peak_current_components_a_m2=abs(j).max(axis=0).tolist(),peak_velocity_components_m_s=abs(response['velocity']).max(axis=0).tolist(),
                minimum_specific_entropy=np.min(u.fluid.entropy,axis=(0,1)).tolist(),maximum_transport_outflow_s=response['maximum_transport_outflow_s'],
                particle_mean_temperature_ev=(np.sum(fluid.particles[...,None]*f['temperature'],axis=(0,1))/fluid.particles.sum()).tolist(),
                minimum_vacuum_map_jacobian=float(response['vacuum_geometry']['determinant'].min()),
                residuals=response['residuals'],pair_power_residual_w=response['pair_power_residual'],weak_charge_residual=response['weak_charge_residual'],
                exchange_coulomb_log_range=response['exchange_coulomb_log_range'],boundary_coordinates=u.fluid.coordinates[domain.boundary].tolist(),
                split_diagnostics=getattr(model,'split_diagnostics',None))
        def arrays(u):return dict(magnetic=u.magnetic,**{k:getattr(u.fluid,k) for k in u.fluid.__dataclass_fields__})
        np.savez_compressed(out/'initial.npz',**arrays(state),particles=fluid.particles,ids=fluid.ids)
        record['initial']=observe(state,0.,initial);save();stride=max(1,steps//80)
        for step in range(1,steps+1):
            state_next=model.step(state,dt,(step-1)*dt);response=model.response(state_next,time=step*dt)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
            external_work+=.5*dt*(previous['external_power']+response['external_power'])
            if integrator=='rk4':source_work+=.5*dt*(previous['source_response_power']+response['source_response_power'])
            heat+=.5*dt*(previous['joule_power']+response['joule_power']);state=state_next;previous=response;record['completed_steps']=step
            raw=out/f'step_{step:05}.npz';np.savez_compressed(raw,**arrays(state))
            if step%stride==0 or step in (1,steps):
                row=observe(state,step*dt,response);row.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(row);save()
            if step==1 or step%(steps//8 or 1)==0:print('material finite coil',nr,step,'/',steps,flush=True)
        record['success']=True
    except Exception as exc:record.update(failure_type=type(exc).__name__,failure_message=str(exc),failed_step=record['completed_steps']+1)
    record['terminal']=True
    if provenance(['examples.material_coil_pulse'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material finite coil terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--steps',type=int,default=2000);p.add_argument('--peak-field',type=float,default=1e-5)
    p.add_argument('--integrator',choices=('rk4','split'),default='rk4')
    a=p.parse_args();d=run(a.destination,nr=a.nr,steps=a.steps,peak_field=a.peak_field,integrator=a.integrator);raise SystemExit(0 if d['success'] else 1)
