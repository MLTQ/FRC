"""Full driven plasma/vacuum evolution against independent cylindrical limits."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, mu_0 as MU0
from plasma_bench.axial_flux_drive import AxialFluxPulse
from plasma_bench.vacuum_cylinder import VacuumCylinderGeometry, VacuumMagnetic
from plasma_bench.vacuum_plasma import VacuumPlasma
from plasma_bench.vacuum_gauss import VacuumGauss
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_plasma import RegularPlasma
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.strong_wall_plasma import StrongWallFluid
from .current_transport_sampling import TransportSampler, bases
from .vacuum_pulse_reference import VacuumPulseReference, VacuumDecayMode
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, *, nr=8, nz=16, steps=160, transport='fixed', seed_amplitude=0.):
    destination=Path(destination)
    if destination.exists() or nr not in (8,16,32) or nz not in (8,16,32) or type(steps) is not int or steps<2 or steps%2:
        raise ValueError('Fresh output, nr/nz8/16/32 and positive even step count required')
    if not np.isfinite(seed_amplitude) or seed_amplitude<0:raise ValueError('Nonnegative finite seed amplitude required')
    started=time.monotonic(); identity=provenance(['examples.vacuum_cylinder_pulse'])
    destination.mkdir(parents=True); freeze(destination,identity)
    radius=.035; length=.24; source_radius=.075; n0=1e21; mi=39.948*atomic_mass; eta0=3.8318200089006576e-5
    drive=AxialFluxPulse(1e-4,12e-6); duration=2*drive.duration_s; dt=duration/steps
    fluid=StrongWallFluid(radius,length,nr//2,nz,mi+ME,degree=5,gammas=(5/3,5/3))
    mesh=MHDMesh(radius,length/2,nr,nz); original_geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    original=RegularPlasma(RegularFieldMap(fluid,original_geometry),mi)
    geometry=VacuumCylinderGeometry(StrongWallGeometry(original_geometry),source_radius=source_radius)
    mapping=RegularFieldMap(fluid,geometry)
    system=VacuumPlasma(mapping,mi,drive=drive,transport=transport,resistivity=eta0); system.original_system=original
    state=np.zeros(fluid.shape); state[...,0]=n0
    state[...,4:]=n0*(np.log([5.,1.])-(fluid.thermal.gammas-1)*np.log(n0/fluid.thermal.n0))
    initial_state=state.copy(); sample=TransportSampler(fluid,mesh,20)
    volume=float(np.sum(sample.weight)); energy_scale=volume*drive.peak_field_t**2/(2*MU0)
    norm=lambda a:float(np.sqrt(np.sum(sample.weight*a*a)/volume))
    options=dict(radius=radius,density=n0,ion_mass=mi,resistivity=eta0,amplitude=drive.peak_field_t,pulse_duration=drive.duration_s)
    reference=VacuumPulseReference(sample.R[:,0],**options,modes=512)
    finer=VacuumPulseReference(sample.R[:,0],**options,modes=1024)
    mode=VacuumDecayMode(reference,length=length,amplitude=seed_amplitude)
    psi=mode.flux(*mesh.corners); psi[:,-1]=psi[:,0]
    field=VacuumMagnetic(mesh,psi,np.zeros((nr,nz)),0.)
    solver=VacuumGauss(system,state,field)
    _,wall_dx,_,_,_=bases(.5*mesh.rf**2,np.array([.5*radius**2]),5)
    exterior_r=np.array([radius,1.5*radius,2*radius]); exterior_z=sample.z
    exterior_maps=np.array([geometry.exterior_samples(np.eye(nz)[j],exterior_r,exterior_z) for j in range(nz)])
    truncation=[]
    for t in (.25*drive.duration_s,.5*drive.duration_s,drive.duration_s,duration):
        a=reference.evaluate(t); b=finer.evaluate(t)
        truncation.append(dict(time=t,bz_max_difference_t=float(np.max(abs(a[0]-b[0]))),j_max_difference_a_m2=float(np.max(abs(a[1]-b[1])))))
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,
        spec=dict(elements=nr//2,points=nz,field_nr=nr,field_nz=nz),transport=transport,drive=asdict(drive),
        radius_m=radius,length_m=length,source_radius_m=source_radius,density_m3=n0,ion_mass_kg=mi,
        gammas=[5/3,5/3],initial_temperature_ev=[5.,1.],fixed_reference_resistivity_ohm_m=eta0,
        seed_amplitude_t=seed_amplitude,seed_radial_wavenumber=mode.q,seed_decay_rate_s=mode.rate,
        pulse_plasma_field_energy_scale_j=energy_scale,source_sha256=identity['source_sha256'],records=[],
        reference_truncation_comparison=truncation,
        scope='Full moving finite-mass plasma in a fixed periodic cylinder, coupled to a magnetic-only exterior extending to infinity and prescribed cylindrical source current at C. No exterior matter, conductivity or free-interface motion. Continuum fixed-eta linear pulse/axial-decay reference; spitzer differences are model differences.')

    def observe(u,b,when):
        model=RegularCurrent(mapping,u[...,0],mi); packed=b.packed_increment(geometry); current=model.current(packed)
        magnetic=sample.magnetic(packed); magnetic[...,2]+=b.axial_bias
        er,et,ez=geometry.split_electric(current)
        old=np.r_[er,et[:-nz],ez]; j=sample.electric(old)
        j[...,1]+=(sample.nx[:,-1]*radius/sample.r)[:,None]*(et[-nz:]@sample.nz.T)[None,:]
        n,v,_=sample.fluid(u)
        entropy=np.stack([sample.scalar(u[...,k])/n for k in (4,5)],axis=-1)
        temperature=np.exp(entropy+(2/3)*np.log(n[...,None]/fluid.thermal.n0))
        energies=dict(thermal=1.5*n*E*np.sum(temperature,axis=-1),kinetic=.5*fluid.mass*n*np.sum(v*v,axis=-1),
            plasma_magnetic=np.sum(magnetic*magnetic,axis=-1)/(2*MU0),current=model.reduced_mass*np.sum(j*j,axis=-1)/(2*n*E**2))
        physical={k:float(np.sum(sample.weight*a)) for k,a in energies.items()}
        components=b.energy_components(model)
        physical['vacuum_induced_magnetic']=components['vacuum_increment']
        physical['vacuum_source_background']=np.pi*(source_radius**2-radius**2)*length*b.axial_bias**2/(2*MU0)
        native=float(np.sum(fluid.volume*fluid.fields(u)['energy'])+sum(components.values())+np.sum(mesh.volumes*model.current_energy(current)))
        expected_b,expected_j=(a[:,None]+np.zeros_like(sample.Z) for a in reference.evaluate(when))
        br,bz,jt=mode.plasma(sample.R,sample.Z,when); expected_b+=bz; expected_j+=jt
        diagnostic=system.derivative(u,b,when).diagnostics
        actual_vacuum=np.tensordot(b.psi_increment[-1,:-1],exterior_maps,axes=(0,0))
        vbr,vbz=mode.vacuum(exterior_r[:,None],exterior_z[None,:],when)
        vac_error=float(np.sqrt(np.mean((actual_vacuum[...,0]-vbr)**2+(actual_vacuum[...,2]-vbz)**2)))
        wall_b=(wall_dx@packed[:(nr+1)*nz].reshape(nr+1,nz))[0]+b.axial_bias
        vacuum_wall=geometry.exterior_samples(b.psi_increment[-1,:-1],np.array([radius]),mesh.zf[:-1])[0,:,2]+b.axial_bias
        return dict(time=when,source_field_t=b.axial_bias,
            mean_plasma_bz_t=float(np.sum(sample.weight*magnetic[...,2])/volume),
            mean_induced_wall_flux=float(np.mean(b.psi_increment[-1,:-1])),
            bz_difference_from_linear_reference_over_drive=norm(magnetic[...,2]-expected_b)/drive.peak_field_t,
            br_difference_from_linear_reference_over_drive=norm(magnetic[...,0]-br)/drive.peak_field_t,
            jtheta_difference_from_linear_reference_over_scale=norm(j[...,1]-expected_j)/(drive.peak_field_t/(MU0*radius)),
            vacuum_vector_difference_from_linear_reference_over_drive=vac_error/drive.peak_field_t,
            reconstructed_interface_bz_jump_over_drive=float(np.sqrt(np.mean((wall_b-vacuum_wall)**2))/drive.peak_field_t),
            bz_min_t=float(magnetic[...,2].min()),bz_max_t=float(magnetic[...,2].max()),maximum_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),
            maximum_speed_m_s=float(np.linalg.norm(v,axis=-1).max()),density_relative_rms_change=norm((n-n0)/n0),
            temperature_mean_ev=(np.sum(sample.weight[...,None]*temperature,axis=(0,1))/volume).tolist(),
            physical_energy_j=physical,native_magnetic_components_j=components,native_total_energy_j=native,
            particle_relative_change=float(np.sum(sample.weight*sample.scalar(u[...,0]-initial_state[...,0]))/(n0*volume)),
            wall_loop_voltage_v=diagnostic['wall_loop_voltage_v'],source_loop_voltage_v=diagnostic['source_loop_voltage_v'],
            source_power_w=diagnostic['source_power_w'],ohmic_power_w=diagnostic['ohmic_power'],
            instantaneous_power_balance_residual_w=diagnostic['power_balance_residual_w'],
            transport_diagnostics={k:diagnostic[k] for k in ('resistivity_ohm_m','resistivity_range_ohm_m','maximum_electron_hall_parameter',
                'maximum_drift_over_species_thermal_speed','weak_field_review','relative_drift_review') if k in diagnostic})

    def save():
        record.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,elapsed_seconds=time.monotonic()-started)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:
            raise RuntimeError('Vacuum pulse exceeded own4GiB/256MiB resource limits')
    np.savez_compressed(destination/'initial.npz',state=state,magnetic=field.packed_increment(geometry),source_field=0.)
    record['initial_sha256']=sha(destination/'initial.npz'); initial=observe(state,field,0.)
    record['initial']=initial; resources(); save(); when=0.; source_work=0.
    for step in range(1,steps+1):
        try:
            updated,evolved,diagnostic=solver.step(state,field,dt,when=when); resources()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error));break
        raw=destination/f'step_{step:04}.npz'
        np.savez_compressed(raw,state=updated,magnetic=evolved.packed_increment(geometry),source_field=evolved.axial_bias,stage_changes=diagnostic.pop('stage_changes'))
        source_work+=diagnostic['source_work_j']; when+=dt; state,field=updated,evolved
        measured=observe(state,field,when)
        measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,integrated_source_work_j=source_work,
            physical_energy_minus_source_over_scale=(sum(measured['physical_energy_j'].values())-sum(initial['physical_energy_j'].values())-source_work)/energy_scale,
            native_energy_minus_source_over_scale=(measured['native_total_energy_j']-initial['native_total_energy_j']-source_work)/energy_scale)
        record['records'].append(measured);record['completed_steps']=step;save()
        if step==1 or step%10==0:print('vacuum pulse',transport,nr,step,'/',steps,'B/reference',measured['bz_difference_from_linear_reference_over_drive'],flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=when)
    if provenance(['examples.vacuum_cylinder_pulse']) != identity:raise RuntimeError('Source changed during vacuum pulse')
    save(); print('vacuum pulse terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--nr',type=int,default=8)
    parser.add_argument('--nz',type=int,default=16);parser.add_argument('--steps',type=int,default=160);parser.add_argument('--transport',choices=('fixed','spitzer'),default='fixed')
    parser.add_argument('--seed-amplitude',type=float,default=0.)
    args=parser.parse_args(); result=run(args.destination,nr=args.nr,nz=args.nz,steps=args.steps,transport=args.transport,seed_amplitude=args.seed_amplitude)
    raise SystemExit(0 if result['success'] else 1)
