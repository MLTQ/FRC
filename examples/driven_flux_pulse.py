"""Full plasma response to one prescribed axial-flux pulse and relaxation."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, mu_0 as MU0
from plasma_bench.axial_flux_drive import AxialFluxPulse, FluxDrivenResistivePlasma, FluxDrivenSpitzerPlasma
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.driven_modal_gauss import DrivenModalGauss
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_plasma import RegularPlasma
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.strong_wall_plasma import StrongWallFluid
from .current_transport_sampling import TransportSampler
from .flux_diffusion_reference import FluxDiffusionReference
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, *, nr=8, steps=160, transport='fixed'):
    destination = Path(destination)
    if destination.exists() or nr not in (8,16,32) or type(steps) is not int or steps < 2 or steps%2:
        raise ValueError('Fresh output, nr8/16/32 and positive even step count required')
    if transport not in ('fixed','spitzer'):
        raise ValueError('Explicit fixed or spitzer transport required')
    started=time.monotonic();identity=provenance(['examples.driven_flux_pulse'])
    destination.mkdir(parents=True);freeze(destination,identity)
    radius=.035;length=.24;n0=1e21;mi=39.948*atomic_mass;eta0=3.8318200089006576e-5
    drive=AxialFluxPulse(1e-4,12e-6);duration=2*drive.duration_s;dt=duration/steps
    fluid=StrongWallFluid(radius,length,nr//2,16,mi+ME,degree=5,gammas=(5/3,5/3))
    mesh=MHDMesh(radius,length/2,nr,16);original_geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    original=RegularPlasma(RegularFieldMap(fluid,original_geometry),mi);geometry=StrongWallGeometry(original_geometry)
    mapping=RegularFieldMap(fluid,geometry)
    system=(FluxDrivenResistivePlasma(mapping,mi,drive=drive,resistivity=eta0,electron_heat_fraction=1.)
        if transport=='fixed' else FluxDrivenSpitzerPlasma(mapping,mi,drive=drive))
    system.original_system=original
    state=np.zeros(fluid.shape);state[...,0]=n0
    state[...,4:]=n0*(np.log([5.,1.])-(fluid.thermal.gammas-1)*np.log(n0/fluid.thermal.n0))
    initial_state=state.copy();field=HarmonicMagnetic(mesh,np.zeros((nr+1,17)),np.zeros((nr,16)),0.)
    solver=DrivenModalGauss(system,state,field);sample=TransportSampler(fluid,mesh,20)
    volume=float(np.sum(sample.weight));energy_scale=volume*drive.peak_field_t**2/(2*MU0)
    norm=lambda a:float(np.sqrt(np.sum(sample.weight*a*a)/volume))
    reference_options=dict(radius=radius,density=n0,ion_mass=mi,resistivity=eta0,
        amplitude=drive.peak_field_t,pulse_duration=drive.duration_s)
    reference=FluxDiffusionReference(sample.R[:,0],**reference_options,modes=512)
    reference_fine=FluxDiffusionReference(sample.R[:,0],**reference_options,modes=1024)
    reference_difference=[]
    for t in (.25*drive.duration_s,.5*drive.duration_s,.75*drive.duration_s,drive.duration_s,duration):
        coarse=reference.evaluate(t);fine=reference_fine.evaluate(t)
        reference_difference.append(dict(time=t,bz_max_difference_t=float(np.max(abs(coarse[0]-fine[0]))),
            jtheta_max_difference_a_m2=float(np.max(abs(coarse[1]-fine[1])))))
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,
        spec=dict(elements=nr//2,points=16,field_nr=nr,field_nz=16),transport=transport,
        drive=asdict(drive),radius_m=radius,length_m=length,density_m3=n0,ion_mass_kg=mi,
        gammas=[5/3,5/3],initial_temperature_ev=[5.,1.],fixed_reference_resistivity_ohm_m=eta0,
        pulse_field_energy_scale_j=energy_scale,source_sha256=identity['source_sha256'],records=[],
        reference_truncation_comparison=reference_difference,
        scope='Full moving plasma driven by enclosed axial flux, not prescribed wall B or finite coils. Fixed-eta continuum reference includes physical current inertia. In spitzer mode reference differences are physical/control differences, not numerical errors. Original current basis has zero wall trace; continuum current need not.')

    def observe(u,b,when):
        model=RegularCurrent(mapping,u[...,0],mi);packed=b.packed_increment(geometry)
        current=model.current(packed);magnetic=sample.magnetic(packed);magnetic[...,2]+=b.axial_bias
        j=sample.electric(current);n,v,_=sample.fluid(u)
        entropy=np.stack([sample.scalar(u[...,k])/n for k in (4,5)],axis=-1)
        temperature=np.exp(entropy+(2/3)*np.log(n[...,None]/fluid.thermal.n0))
        energies=dict(thermal=1.5*n*E*np.sum(temperature,axis=-1),kinetic=.5*fluid.mass*n*np.sum(v*v,axis=-1),
            magnetic=np.sum(magnetic*magnetic,axis=-1)/(2*MU0),
            current=model.reduced_mass*np.sum(j*j,axis=-1)/(2*n*E**2))
        physical={k:float(np.sum(sample.weight*a)) for k,a in energies.items()}
        native=float(np.sum(fluid.volume*fluid.fields(u)['energy'])+.5*packed@(model.hessian@packed)
            +np.sum(mesh.volumes*model.current_energy(current))+volume*b.axial_bias**2/(2*MU0))
        expected_b,expected_j=(a[:,None]+np.zeros_like(sample.Z) for a in reference.evaluate(when))
        diagnostic=system.derivative(u,b,when).diagnostics
        return dict(time=when,mean_flux_field_t=b.axial_bias,
            bz_difference_from_fixed_eta_reference_over_drive=norm(magnetic[...,2]-expected_b)/drive.peak_field_t,
            jtheta_difference_from_fixed_eta_reference_over_scale=norm(j[...,1]-expected_j)/(drive.peak_field_t/(MU0*radius)),
            bz_min_t=float(magnetic[...,2].min()),bz_max_t=float(magnetic[...,2].max()),
            maximum_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),
            maximum_speed_m_s=float(np.linalg.norm(v,axis=-1).max()),
            maximum_radial_velocity_m_s=float(abs(v[...,0]).max()),
            minimum_density_m3=float(n.min()),density_relative_rms_change=norm((n-n0)/n0),
            temperature_mean_ev=(np.sum(sample.weight[...,None]*temperature,axis=(0,1))/volume).tolist(),
            physical_energy_j=physical,native_total_energy_j=native,
            particle_relative_change=float(np.sum(sample.weight*sample.scalar(u[...,0]-initial_state[...,0]))/(n0*volume)),
            physical_mean_B_minus_prescribed_t=float(np.sum(sample.weight*magnetic[...,2])/volume-b.axial_bias),
            wall_loop_voltage_v=diagnostic['wall_loop_voltage_v'],source_power_w=diagnostic['source_power_w'],
            source_plasma_power_w=diagnostic['source_plasma_power_w'],ohmic_power_w=diagnostic['ohmic_power'],
            instantaneous_power_balance_residual_w=diagnostic['power_balance_residual_w'],
            transport_diagnostics={k:diagnostic[k] for k in ('resistivity_ohm_m','resistivity_range_ohm_m',
                'maximum_electron_hall_parameter','maximum_drift_over_species_thermal_speed',
                'weak_field_review','relative_drift_review') if k in diagnostic})

    def save():
        record.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.monotonic()-started)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')

    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 4*1024**3 or shutil.disk_usage('.').free < 256*1024**2:
            raise RuntimeError('Driven pulse exceeded own4GiB/256MiB resource limits')

    np.savez_compressed(destination/'initial.npz',state=state,magnetic=field.packed_increment(geometry),mean_flux_field=0.)
    record['initial_sha256']=sha(destination/'initial.npz');initial=observe(state,field,0.)
    record['initial']=initial;resources();save();when=0.;source_work=0.
    for step in range(1,steps+1):
        try:
            updated,evolved,diagnostic=solver.step(state,field,dt,when=when);resources()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error));break
        raw=destination/f'step_{step:04}.npz'
        np.savez_compressed(raw,state=updated,magnetic=evolved.packed_increment(geometry),
            mean_flux_field=evolved.axial_bias,stage_changes=diagnostic.pop('stage_changes'))
        source_work+=diagnostic['source_work_j'];when+=dt;state,field=updated,evolved
        measured=observe(state,field,when)
        measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,
            integrated_source_work_j=source_work,
            physical_energy_minus_source_over_scale=(sum(measured['physical_energy_j'].values())-sum(initial['physical_energy_j'].values())-source_work)/energy_scale,
            native_energy_minus_source_over_scale=(measured['native_total_energy_j']-initial['native_total_energy_j']-source_work)/energy_scale)
        record['records'].append(measured);record['completed_steps']=step;save()
        if step==1 or step%10==0:
            print('flux pulse',transport,nr,step,'/',steps,'B/reference',measured['bz_difference_from_fixed_eta_reference_over_drive'],flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=when)
    if provenance(['examples.driven_flux_pulse']) != identity:
        raise RuntimeError('Source changed during driven pulse')
    save();print('flux pulse terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--nr',type=int,default=8)
    parser.add_argument('--steps',type=int,default=160);parser.add_argument('--transport',choices=('fixed','spitzer'),default='fixed')
    args=parser.parse_args();run(args.destination,nr=args.nr,steps=args.steps,transport=args.transport)
