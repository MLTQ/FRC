"""Complete full-vector material plasma expansion with independent force/virial comparisons."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,electron_mass as ME,mu_0 as MU0
from scipy.sparse.linalg import splu
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_vector_space import MaterialVectorSpace
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_vector_current import MaterialVectorCurrent
from .implicit_plasma_records import provenance,freeze,sha


def compact_vector(r,z,radius,length,bp,bt,density):
    """Analytic compact psi=(Bp/2)*r²*hr^4*hz^4 and Btheta=Bt*r/R*hr³*hz³."""
    s=(r/radius)**2;t=(z/length)**2;hr=1-s;hz=1-t
    b=np.stack((4*bp*r*z/length**2*hr**4*hz**3,bt*r/radius*hr**3*hz**3,bp*hr**3*hz**4*(1-5*s)),axis=-1)
    j=np.stack((6*bt*r*z/(MU0*radius*length**2)*hr**3*hz**2,
        bp*r/MU0*(4*hr**4*hz**2*(1-7*t)/length**2+8*hr**2*(2-5*s)*hz**4/radius**2),
        bt/(MU0*radius)*hz**3*(2*hr**3-6*s*hr**2)),axis=-1)
    return b,j,density*(hr*hz)**1.5


def initial_force_comparison(model,state,radius,length,bp,bt,density):
    f=model.fluid;r=f.reference_geometry['points'][...,0].ravel();z=f.reference_geometry['points'][...,1].ravel()
    b,j,n=compact_vector(r,z,radius,length,bp,bt,density)
    def stress(rr,zz):
        _,jj,nn=compact_vector(rr,zz,radius,length,bp,bt,density)
        return model.mu/E**2*jj[..., :,None]*jj[...,None,:]/nn[...,None,None]
    epsilon=1e-25;radial=stress(r+1j*epsilon,z).imag/epsilon;axial=stress(r,z+1j*epsilon).imag/epsilon;tensor=stress(r,z)
    divergence=np.column_stack((radial[:,0,0]+axial[:,0,2]+(tensor[:,0,0]-tensor[:,1,1])/r,
        radial[:,1,0]+axial[:,1,2]+2*tensor[:,1,0]/r,radial[:,2,0]+axial[:,2,2]+tensor[:,2,0]/r))
    hr=1-(r/radius)**2;hz=1-(z/length)**2
    pressure=10*E*n[:,None]*np.column_stack((r*hz/radius**2,np.zeros_like(r),z*hr/length**2))
    exact=(pressure+np.cross(j,b)-divergence)/(f.mass*n[:,None])
    rhs=model.derivative(state);acceleration=np.zeros_like(state.fluid.coordinates)
    acceleration[f.active_r,0]=f.mass_factors[0].solve(rhs.fluid.momentum[f.active_r,0]);acceleration[:,1]=f.mass_factors[1].solve(rhs.fluid.momentum[:,1])
    poloidal=np.einsum('qn,cna->cqa',f.shape,acceleration[f.ids]).reshape(-1,2)
    angular_factor=splu(f.gram(f.mass*f.particles*f.reference_geometry['points'][...,0]**2))
    omega_dot=angular_factor.solve(rhs.fluid.angular_momentum)
    toroidal=r*np.einsum('qn,cn->cq',f.shape,omega_dot[f.ids]).ravel()
    actual=np.column_stack((poloidal[:,0],toroidal,poloidal[:,1]));weights=f.particles.ravel();scale=10*E/(f.mass*radius)
    return dict(acceleration_particle_rms_error_over_C_R=float(np.sqrt(np.sum(weights*np.sum((actual-exact)**2,axis=-1))/weights.sum())/scale),
        component_acceleration_rms_error_over_C_R=(np.sqrt(np.sum(weights[:,None]*(actual-exact)**2,axis=0)/weights.sum())/scale).tolist())


def run(destination,*,nr=4,nz=8,degree=3,quadrature=6,steps=2000,resistivity=.001,bp=2e-4,bt=5e-5):
    out=Path(destination)
    if out.exists() or degree not in (1,3,5) or min(nr,nz)<degree or quadrature<2*degree or steps<4:
        raise ValueError('Fresh complete vector trajectory and supported resolution required')
    started=time.monotonic();identity=provenance(['examples.material_vector_expansion']);out.mkdir(parents=True);freeze(out,identity)
    radius=.035;length=.07;mi=39.948*atomic_mass;nc=1e18;mesh=MHDMesh(radius,length,nr,nz)
    domain=MaterialCylinder(mesh,quadrature=quadrature);space=MaterialVectorSpace(mesh,degree=degree,quadrature=quadrature)
    r=space.reference_r.reshape(nr*nz,-1);z=space.reference_z.reshape(nr*nz,-1);profile=(1-(r/radius)**2)*(1-(z/length)**2)
    b,j,n=compact_vector(r,z,radius,length,bp,bt,nc)
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,n,np.broadcast_to(profile[...,None],profile.shape+(2,)),mi)
    fluid.initial.entropy[:]=-(2/3)*np.log(nc);model=MaterialVectorCurrent(fluid,space,resistivity=resistivity)
    state=model.project_initial(fluid.initial,b,j);initial=model.response(state);duration=radius/np.sqrt(10*E/fluid.mass);dt=duration/steps
    energy0=initial['total_energy'];inertia0=float(np.sum(fluid.mass*fluid.particles*np.sum(initial['fields']['points']**2,axis=-1)))
    volume=initial['volume'];current=E*initial['fields']['n'].ravel()[:,None]*initial['drift']
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,nr=nr,nz=nz,degree=degree,quadrature=quadrature,
        radius_m=radius,half_length_m=length,ion_mass_kg=mi,central_density_m3=nc,poloidal_amplitude_t=bp,toroidal_amplitude_t=bt,
        resistivity_ohm_m=resistivity,coulomb_exchange=True,entropy_transport='upwind',source_sha256=identity['source_sha256'],initial_total_energy_j=energy0,records=[],
        initial_magnetic_relative_l2_error=float(np.sqrt(np.sum(volume[:,None]*(initial['magnetic_field']-b.reshape(-1,3))**2)/np.sum(volume[:,None]*b.reshape(-1,3)**2))),
        initial_current_relative_l2_error=float(np.sqrt(np.sum(volume[:,None]*(current-j.reshape(-1,3))**2)/np.sum(volume[:,None]*j.reshape(-1,3)**2))),
        initial_force_comparison=initial_force_comparison(model,state,radius,length,bp,bt,nc),
        scope='All axisymmetric B/J/flow components, actual material rotation, full finite-mass Hall/current stress and species heat. Compact B=0 boundary; exact empty exterior. Initialized flux is not coil-driven FRC formation; no physical edge/transport validity claim.')
    wb_integral=0.;twb_integral=0.;heat_integral=0.;previous=initial
    def observe(u,when,a):
        fields=a['fields'];n=fields['n'].ravel();j=E*n[:,None]*a['drift'];b=a['magnetic_field']
        inertia=float(np.sum(fluid.mass*fluid.particles*np.sum(fields['points']**2,axis=-1)))
        predicted=inertia0+2*energy0*when**2-2*(when*wb_integral-twb_integral)
        return dict(time=when,total_energy_j=a['total_energy'],relative_total_energy_change=a['total_energy']/energy0-1,
            magnetic_energy_j=a['magnetic_energy'],current_energy_j=a['current_energy'],bulk_kinetic_energy_j=fields['kinetic_energy'],internal_energy_j=fields['internal_energy'].tolist(),
            joule_power_w=a['joule_power'],joule_heat_j=heat_integral,peak_b_components_t=abs(b).max(axis=0).tolist(),peak_current_components_a_m2=abs(j).max(axis=0).tolist(),
            peak_b_t=float(np.linalg.norm(b,axis=-1).max()),peak_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),peak_drift_m_s=float(np.linalg.norm(a['drift'],axis=-1).max()),
            peak_velocity_components_m_s=abs(a['velocity']).max(axis=0).tolist(),peak_material_angle_rad=float(abs(u.fluid.angle).max()),
            total_angular_momentum_kg_m2_s=float(u.fluid.angular_momentum.sum()),total_axial_momentum_kg_m_s=float(u.fluid.momentum[:,1].sum()),
            minimum_specific_entropy=np.min(u.fluid.entropy,axis=(0,1)).tolist(),maximum_transport_outflow_s=a['maximum_transport_outflow_s'],
            minimum_density_m3=float(n.min()),minimum_temperature_ev=float(fields['temperature'].min()),particle_mean_temperature_ev=(np.sum(fluid.particles[...,None]*fields['temperature'],axis=(0,1))/fluid.particles.sum()).tolist(),
            virial_inertia_error_over_initial=(inertia-predicted)/inertia0,residuals=a['residuals'],weak_charge_residual=a['weak_charge_residual'],pair_power_residual_w=a['pair_power_residual'],
            exchange_coulomb_log_range=a['exchange_coulomb_log_range'],boundary_coordinates=u.fluid.coordinates[domain.boundary].tolist())
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
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
        except Exception as exc:
            record.update(failed_step=step,failure_type=type(exc).__name__,failure_message=str(exc));break
        state=updated;when=step*dt;wb_integral+=.5*dt*(previous['magnetic_energy']+response['magnetic_energy'])
        twb_integral+=.5*dt*(t0*previous['magnetic_energy']+when*response['magnetic_energy']);heat_integral+=.5*dt*(previous['joule_power']+response['joule_power']);previous=response
        raw=out/f'step_{step:05}.npz';np.savez_compressed(raw,**arrays(state));record['completed_steps']=step
        if step%stride==0 or step==1 or step==steps:
            row=observe(state,when,response);row.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw));record['records'].append(row);save()
        if step==1 or step%(steps//8 or 1)==0:print('material vector',nr,step,'/',steps,flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.material_vector_expansion'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material vector terminal',record['success'],record.get('failure_message',''),flush=True);return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--nz',type=int,default=8)
    p.add_argument('--degree',type=int,default=3);p.add_argument('--quadrature',type=int,default=6);p.add_argument('--steps',type=int,default=2000)
    p.add_argument('--resistivity',type=float,default=.001);p.add_argument('--bp',type=float,default=2e-4);p.add_argument('--bt',type=float,default=5e-5)
    a=p.parse_args();d=run(a.destination,nr=a.nr,nz=a.nz,degree=a.degree,quadrature=a.quadrature,steps=a.steps,resistivity=a.resistivity,bp=a.bp,bt=a.bt)
    raise SystemExit(0 if d['success'] else 1)
