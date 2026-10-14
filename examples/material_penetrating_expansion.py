"""Simultaneous heated plasma/current expansion with a penetrating magnetic vacuum."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_magnetic_vacuum import MaterialMagneticVacuum
from plasma_bench.material_vacuum_space import MaterialVacuumSpace
from plasma_bench.material_vacuum_current import MaterialVacuumCurrent
from .material_vacuum_motion import reference as exterior_reference
from .implicit_plasma_records import provenance,freeze,sha


def spherical_current(r,z,radius,central_field):
    """Smooth compact volume current, analytic dipole beyond its support."""
    square=r*r+z*z;u=square/radius**2;moment=central_field*radius**3/8.75
    inside=u<=1;f=np.empty_like(r);df=np.empty_like(r);j=np.zeros(r.shape+(3,))
    f[inside]=moment/radius**3*(35-42*u[inside]+15*u[inside]**2)/8
    df[inside]=moment/radius**5*(-42+30*u[inside])/8
    f[~inside]=moment/square[~inside]**1.5;df[~inside]=-1.5*moment/square[~inside]**2.5
    j[...,1][inside]=52.5*moment/(MU0*radius**5)*r[inside]*(1-u[inside])
    b=np.stack((-2*r*z*df,np.zeros_like(r),2*f+2*r*r*df),axis=-1)
    return r*r*f,b,j


def force_comparison(model,state,b,j,n,r,z,R,L):
    f=model.fluid;hr=1-(r/R)**2;hz=1-(z/L)**2
    force=10*E*n[:,None]*np.column_stack((r*hz/R**2,np.zeros_like(r),z*hr/L**2))+np.cross(j,b)
    # Pure azimuthal current: div(mu*n*w*w)_r=-mu*n*wtheta²/r.
    force[:,0]+=model.mu*j[:,1]**2/(E**2*n*r)
    target=force[:,[0,2]]/(f.mass*n[:,None]);rhs=model.derivative(state)
    nodes=np.zeros_like(state.fluid.coordinates)
    nodes[f.active_r,0]=f.mass_factors[0].solve(rhs.fluid.momentum[f.active_r,0]);nodes[:,1]=f.mass_factors[1].solve(rhs.fluid.momentum[:,1])
    actual=np.einsum('qn,cna->cqa',f.shape,nodes[f.ids]).reshape(-1,2)
    return float(np.sqrt(np.sum(f.particles.ravel()*np.sum((actual-target)**2,axis=1))/f.particles.sum())/(10*E/(f.mass*R)))


def run(destination,*,nr=4,degree=3,quadrature=6,vacuum_cells=4,modes=6,steps=2000,resistivity=.001):
    out=Path(destination)
    if out.exists() or steps<4 or quadrature<2*degree:raise ValueError('Fresh complete trajectory and resolved field quadrature required')
    started=time.monotonic();identity=provenance(['examples.material_penetrating_expansion']);out.mkdir(parents=True);freeze(out,identity)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,nr=nr,nz=2*nr,degree=degree,quadrature=quadrature,vacuum_cells=vacuum_cells,modes=modes,
        source_sha256=identity['source_sha256'],radius_m=.035,half_length_m=.07,central_density_m3=1e18,central_field_t=2e-4,
        resistivity_ohm_m=resistivity,entropy_transport='upwind',coulomb_exchange=True,records=[],
        scope='Full moving finite-mass plasma/current/species and penetrating magnetic vacuum through both ends and side. Initialized compact spherical volume current, prescribed resistance; no coils, formation, physical heated-edge or local-current convergence claim.')
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        R=record['radius_m'];L=record['half_length_m'];nc=record['central_density_m3'];a=.7*R;mesh=MHDMesh(R,L,nr,2*nr)
        domain=MaterialCylinder(mesh,quadrature=quadrature);space=MaterialVacuumSpace(mesh,degree=degree,quadrature=quadrature)
        vacuum=MaterialMagneticVacuum(mesh,degree=degree,quadrature_order=quadrature,vacuum_cells=vacuum_cells,modes=modes)
        r=space.reference_r;z=space.reference_z;profile=(1-(r/R)**2)*(1-(z/L)**2);density=nc*profile**1.5
        shape=(nr*2*nr,-1);n=density.reshape(shape);temperature=np.repeat(profile.reshape(shape)[...,None],2,axis=-1)
        fluid=MaterialTwoSpecies(domain,domain.fluid_reference,n,temperature,39.948*atomic_mass);fluid.initial.entropy[:]=-(2/3)*np.log(nc)
        model=MaterialVacuumCurrent(fluid,space,vacuum,resistivity=resistivity)
        _,b,j=spherical_current(r,z,a,record['central_field_t']);x=fluid.initial.coordinates[vacuum.boundary_nodes]
        boundary,_,_=spherical_current(x[:,0],x[:,1],a,record['central_field_t'])
        state=model.project_initial(fluid.initial,b,j,boundary);initial=model.response(state)
        duration=R/np.sqrt(10*E/fluid.mass);dt=duration/steps;energy0=initial['total_energy']
        inertia0=float(np.sum(fluid.mass*fluid.particles*np.sum(initial['fields']['points']**2,axis=-1)))
        exact_vac=exterior_reference(fluid.initial.coordinates,np.zeros_like(fluid.initial.coordinates),domain.ids,record['central_field_t']*a**3/8.75)
        weight=initial['volume'];current=E*n.ravel()[:,None]*initial['drift']
        record.update(dt=dt,duration=duration,ion_mass_kg=fluid.ion_mass,current_support_radius_m=a,
            initial_total_energy_j=energy0,initial_magnetic_relative_l2_error=float(np.sqrt(np.sum(weight[:,None]*(initial['magnetic_field']-b)**2)/np.sum(weight[:,None]*b*b))),
            initial_current_relative_l2_error=float(np.sqrt(np.sum(weight[:,None]*(current-j)**2)/np.sum(weight[:,None]*j*j))),
            initial_vacuum_energy_relative_error=initial['vacuum_magnetic_energy']/exact_vac['energy']-1,
            initial_acceleration_rms_error_over_C_R=force_comparison(model,state,b,j,density,r,z,R,L))
        wb_integral=0.;twb_integral=0.;heat=0.;previous=initial
        def observe(u,t,response):
            fields=response['fields'];density=fields['n'].ravel();current=E*density[:,None]*response['drift']
            inertia=float(np.sum(fluid.mass*fluid.particles*np.sum(fields['points']**2,axis=-1)))
            predicted=inertia0+2*energy0*t*t-2*(t*wb_integral-twb_integral)
            return dict(time=t,total_energy_j=response['total_energy'],relative_total_energy_change=response['total_energy']/energy0-1,
                magnetic_energy_j=response['magnetic_energy'],vacuum_magnetic_energy_j=response['vacuum_magnetic_energy'],current_energy_j=response['current_energy'],
                peak_b_components_t=abs(response['magnetic_field']).max(axis=0).tolist(),peak_current_components_a_m2=abs(current).max(axis=0).tolist(),
                peak_velocity_components_m_s=abs(response['velocity']).max(axis=0).tolist(),peak_drift_m_s=float(np.linalg.norm(response['drift'],axis=1).max()),
                internal_energy_j=fields['internal_energy'].tolist(),bulk_kinetic_energy_j=fields['kinetic_energy'],joule_heat_j=heat,
                particle_mean_temperature_ev=(np.sum(fluid.particles[...,None]*fields['temperature'],axis=(0,1))/fluid.particles.sum()).tolist(),
                minimum_specific_entropy=np.min(u.fluid.entropy,axis=(0,1)).tolist(),maximum_transport_outflow_s=response['maximum_transport_outflow_s'],
                minimum_vacuum_map_jacobian=float(response['vacuum_geometry']['determinant'].min()),
                total_axial_momentum_kg_m_s=float(u.fluid.momentum[:,1].sum()),total_angular_momentum_kg_m2_s=float(u.fluid.angular_momentum.sum()),
                virial_inertia_error_over_initial=(inertia-predicted)/inertia0,residuals=response['residuals'],pair_power_residual_w=response['pair_power_residual'],
                weak_charge_residual=response['weak_charge_residual'],exchange_coulomb_log_range=response['exchange_coulomb_log_range'])
        def arrays(u):return dict(magnetic=u.magnetic,**{k:getattr(u.fluid,k) for k in u.fluid.__dataclass_fields__})
        np.savez_compressed(out/'initial.npz',**arrays(state),particles=fluid.particles,ids=fluid.ids)
        record['initial']=observe(state,0.,initial);save();stride=max(1,steps//80)
        for step in range(1,steps+1):
            updated=model.step(state,dt);response=model.response(updated)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
            t=step*dt;wb_integral+=.5*dt*(previous['magnetic_energy']+response['magnetic_energy'])
            twb_integral+=.5*dt*((t-dt)*previous['magnetic_energy']+t*response['magnetic_energy']);heat+=.5*dt*(previous['joule_power']+response['joule_power'])
            state=updated;previous=response;record['completed_steps']=step
            raw=out/f'step_{step:05}.npz';np.savez_compressed(raw,**arrays(state))
            if step%stride==0 or step in (1,steps):
                row=observe(state,t,response);row.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(row);save()
            if step==1 or step%(steps//8 or 1)==0:print('material vacuum plasma',nr,step,'/',steps,flush=True)
        record['success']=True
    except Exception as exc:record.update(failure_type=type(exc).__name__,failure_message=str(exc),failed_step=record['completed_steps']+1)
    record['terminal']=True
    if provenance(['examples.material_penetrating_expansion'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material vacuum plasma terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--degree',type=int,default=3);p.add_argument('--quadrature',type=int,default=6)
    p.add_argument('--vacuum-cells',type=int,default=4);p.add_argument('--modes',type=int,default=6);p.add_argument('--steps',type=int,default=2000);p.add_argument('--resistivity',type=float,default=.001)
    a=p.parse_args();d=run(a.destination,nr=a.nr,degree=a.degree,quadrature=a.quadrature,vacuum_cells=a.vacuum_cells,modes=a.modes,steps=a.steps,resistivity=a.resistivity)
    raise SystemExit(0 if d['success'] else 1)
