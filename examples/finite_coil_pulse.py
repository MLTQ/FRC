"""Full finite-tube coil pulse against the independent weak-conductivity limit."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass,electron_mass as ME,elementary_charge as E,mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_cylinder_fluid import FiniteCylinderFluid
from plasma_bench.finite_cylinder_geometry import FiniteCylinderGeometry
from plasma_bench.finite_cylinder_map import FiniteCylinderMap
from plasma_bench.finite_coil_source import FiniteCoilSource,FiniteCoilMagnetic
from plasma_bench.finite_coil_plasma import FiniteCoilPlasma
from plasma_bench.finite_coil_gauss import FiniteCoilGauss
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.waveforms import Unipolar
from .implicit_plasma_records import provenance,freeze,sha


def line_integral_potential(r,z,coils,order=256):
    """Direct azimuthal Biot-Savart Atheta integral, no elliptic field routine."""
    r,z=np.broadcast_arrays(r,z);rr=r.ravel();zz=z.ravel();result=np.zeros(rr.shape)
    q,w=np.polynomial.legendre.leggauss(order);phi=np.pi*(q+1);weight=np.pi*w;cos=np.cos(phi)
    for start in range(0,len(rr),2048):
        rv=rr[start:start+2048,None];zv=zz[start:start+2048,None]
        for a,zc,turns in coils:
            distance=np.sqrt(rv*rv+a*a-2*a*rv*cos+(zv-zc)**2)
            result[start:start+2048]+=turns*MU0*a/(4*np.pi)*np.sum(weight*cos/distance,axis=1)
    return result.reshape(r.shape)


def run(destination,*,nr=8,nz=16,fluid_z=16,steps=80,resistivity=1.,transport='fixed',peak_field=1e-4):
    destination=Path(destination)
    if destination.exists() or nr<6 or nz<6 or fluid_z<6 or steps<4 or steps%2 or resistivity<=0:
        raise ValueError('Fresh destination, resolved grids, positive resistance and even full-pulse steps required')
    if not np.isfinite(peak_field) or peak_field<=0:raise ValueError('Positive finite prescribed peak field required')
    started=time.monotonic();identity=provenance(['examples.finite_coil_pulse']);destination.mkdir(parents=True);freeze(destination,identity)
    radius=.075;length=1.;density=1e20;mi=39.948*atomic_mass;temperature=np.array([5.,1.]);pulse=12e-6;duration=24e-6;dt=duration/steps
    fluid=FiniteCylinderFluid(radius,length,nr//2,fluid_z,mi+ME)
    mesh=MHDMesh(radius,length/2,nr,nz);g=FiniteCylinderGeometry(mesh,quadrature_order=12,vacuum_cells=12,modes=8)
    mapping=FiniteCylinderMap(fluid,g);coils=((.12,-.15,1.),(.12,.15,1.));source=FiniteCoilSource(g,coils)
    center_unit=sum(MU0*a*a*turns/(2*(a*a+z*z)**1.5) for a,z,turns in coils)
    drive=Unipolar(peak_field/center_unit,rise=pulse/2,hold=0.,fall=pulse/2)
    system=FiniteCoilPlasma(mapping,mi,source=source,drive=drive,transport=transport,resistivity=resistivity)
    state=np.zeros(fluid.shape);state[...,0]=density
    state[...,4:]=density*(np.log(temperature)-(2/3)*np.log(density/fluid.thermal.n0));initial_state=state.copy();f0=fluid.fields(state)
    field=FiniteCoilMagnetic(source,np.zeros((nr+1,nz+1)),np.zeros((nr,nz)),0.)
    solver=FiniteCoilGauss(system,state,field)
    reference_a=line_integral_potential(g.r,g.z,coils);finer_a=line_integral_potential(g.r,g.z,coils,384)
    volume=float(np.sum(g.weight));rms=lambda v:float(np.sqrt(np.sum(g.weight*v*v)/volume))
    current_scale=rms(reference_a)*drive.amplitude*np.pi/pulse/resistivity
    expected_heat=float(np.sum(g.weight*reference_a**2)/resistivity*np.pi**2*drive.amplitude**2/(2*pulse))
    energy_scale=volume*peak_field**2/(2*MU0)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,
        field_nr=nr,field_nz=nz,fluid_axial_points=fluid_z,fluid_radial_elements=nr//2,transport=transport,
        weak_conductivity_reference_applicable=transport=='fixed',fixed_resistivity_ohm_m=resistivity,radius_m=radius,length_m=length,density_m3=density,ion_mass_kg=mi,
        initial_temperature_ev=temperature.tolist(),coils=coils,drive=asdict(drive),peak_center_source_field_t=peak_field,
        reference_current_rms_peak_a_m2=current_scale,weak_conductivity_expected_heat_j=expected_heat,
        field_energy_scale_j=energy_scale,source_sha256=identity['source_sha256'],
        direct_integral_order_difference_relative=rms(reference_a-finer_a)/rms(finer_a),
        analytic_source_A_difference_relative=rms(reference_a-source.a_unit[...,1])/rms(reference_a),records=[],
        reference_scope='Weak fixed-resistance limit only; not a Spitzer reference. Finite induction/inertia contributions remain.',
        scope='Actual moving finite-mass plasma and magnetic vacuum around sealed finite cylinder. Two circular coils with prescribed current. Fixed eta=1 Ohm m is a weak-conductivity mathematical limit, not laboratory transport. Leading Jtheta=-Idot*Aunit/eta; finite magnetic response/inertia and nonlinear changes remain in simulation. Filament self-energy baseline excluded; mutual energy and source response work retained.')
    coil_r=np.array([c[0] for c in coils]);coil_z=np.array([c[1] for c in coils]);coil_turns=np.array([c[2] for c in coils])
    # Linear vacuum maps reused for observations; no extra exterior solve per step.
    vacuum_link=np.array([2*np.pi*np.sum(coil_turns*g.vacuum.samples(np.eye(len(g.boundary_indices))[k],coil_r,coil_z)[0]) for k in range(len(g.boundary_indices))])

    def observe(u,b,t):
        f=fluid.fields(u);model=RegularCurrent(mapping,u[...,0],mi);packed=b.packed_increment(g);current=model.current(packed);j=g.sample_electric(current)
        magnetic=b.samples(g);velocity=mapping.velocity(u);nt=mapping.scalar(u[...,0]);entropy=mapping.scalar(u[...,4])/nt
        temps=np.stack([np.exp(mapping.scalar(u[...,k])/nt+(2/3)*np.log(nt/fluid.thermal.n0)) for k in (4,5)],axis=-1)
        expected=-float(drive.derivative(t))*reference_a/resistivity
        dn=u[...,0]-initial_state[...,0];ratio=dn/initial_state[...,0]
        ds=(u[...,4:]-initial_state[...,4:]-f0['entropy']*dn[...,None])/u[...,0,None]
        dh=ds+(2/3)*np.log1p(ratio)[...,None]
        internal=np.sum(initial_state[...,0,None]*f0['specific_internal']*(ratio[...,None]+(1+ratio[...,None])*np.expm1(dh)),axis=-1)
        thermal=float(np.sum(fluid.volume*internal));kinetic=float(np.sum(fluid.volume*.5*fluid.mass*u[...,0]*np.sum(f['velocity']**2,axis=-1)))
        linkage=source.linkage(current);vacuum_linkage=float(vacuum_link@(g.boundary_map@packed))
        energies=dict(thermal_change=thermal,bulk_kinetic=kinetic,current_kinetic=float(np.sum(mesh.volumes*model.current_energy(current))),
            induced_magnetic=float(.5*packed@(g.hessian@packed)),source_mutual=b.source_amplitude*linkage)
        diagnostic=system.derivative(u,b,t).diagnostics
        return dict(time=t,source_current_a=b.source_amplitude,source_center_bz_t=b.source_amplitude*center_unit,
            current_error_over_reference_peak=rms(j[...,1]-expected)/current_scale if transport=='fixed' else None,
            current_rms_a_m2=rms(j[...,1]),reference_current_rms_a_m2=rms(expected) if transport=='fixed' else None,
            maximum_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),maximum_speed_m_s=float(np.linalg.norm(velocity,axis=-1).max()),
            bz_min_t=float(magnetic[...,2].min()),bz_max_t=float(magnetic[...,2].max()),
            mean_plasma_bz_t=float(np.sum(g.weight*magnetic[...,2])/volume),
            temperature_mean_ev=(np.sum(g.weight[...,None]*temps,axis=(0,1,2,3))/volume).tolist(),
            particle_relative_change=float(np.sum(fluid.volume*dn)/(density*volume)),
            energy_change_components_j=energies,response_energy_change_j=sum(energies.values()),
            source_linked_flux_wb=linkage,vacuum_reconstructed_linkage_wb=vacuum_linkage,
            source_power_w=diagnostic['source_power_w'],ohmic_power_w=diagnostic['ohmic_power'],
            instantaneous_power_balance_residual_w=diagnostic['power_balance_residual_w'],
            induction_residuals=diagnostic['induction_residuals'],
            transport_diagnostics={k:diagnostic[k] for k in ('resistivity_range_ohm_m','maximum_electron_hall_parameter',
                'maximum_drift_over_species_thermal_speed','weak_field_review','relative_drift_review') if k in diagnostic})
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    def guard():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:
            raise RuntimeError('Finite-coil run exceeded own4GiB/256MiB resource bounds')
    np.savez_compressed(destination/'initial.npz',state=state,magnetic=field.packed_increment(g),r=g.r,z=g.z,weight=g.weight,source_a=reference_a)
    record['initial_sha256']=sha(destination/'initial.npz');record['initial']=observe(state,field,0.);save();guard();when=0.;work=0.
    for step in range(1,steps+1):
        try:updated,evolved,diagnostic=solver.step(state,field,dt,when=when);guard()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error));break
        raw=destination/f'step_{step:04}.npz';stage=diagnostic.pop('stage_changes')
        np.savez_compressed(raw,state=updated,magnetic=evolved.packed_increment(g),source_current=evolved.source_amplitude,stage_changes=stage)
        when+=dt;state,field=updated,evolved;work+=diagnostic['source_work_j']
        measured=observe(state,field,when)
        measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,integrated_source_work_j=work,
            energy_minus_source_over_field_scale=(measured['response_energy_change_j']-work)/energy_scale,
            energy_minus_source_over_reference_heat=(measured['response_energy_change_j']-work)/expected_heat if transport=='fixed' else None)
        record['records'].append(measured);record['completed_steps']=step;save()
        if step==1 or step%10==0:print('finite coil',nr,nz,step,'/',steps,'current error',measured['current_error_over_reference_peak'],flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=when)
    if provenance(['examples.finite_coil_pulse'])!=identity:raise RuntimeError('Source changed during finite-coil pulse')
    save();print('finite coil terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--nr',type=int,default=8);parser.add_argument('--nz',type=int,default=16)
    parser.add_argument('--fluid-z',type=int,default=16);parser.add_argument('--steps',type=int,default=80);parser.add_argument('--resistivity',type=float,default=1.)
    parser.add_argument('--peak-field',type=float,default=1e-4)
    parser.add_argument('--transport',choices=('fixed','spitzer'),default='fixed');args=parser.parse_args()
    result=run(args.destination,nr=args.nr,nz=args.nz,fluid_z=args.fluid_z,steps=args.steps,resistivity=args.resistivity,transport=args.transport,peak_field=args.peak_field)
    raise SystemExit(0 if result['success'] else 1)
