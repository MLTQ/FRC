"""Full moving, nonuniform temperature relaxation over one periodic transit."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, mu_0 as MU0
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.coulomb_exchange_wall_plasma import CoulombExchangeWallPlasma
from plasma_bench.guarded_modal_gauss import GuardedModalGauss
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_plasma import RegularPlasma
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.strong_wall_plasma import StrongWallFluid
from .coulomb_temperature_relaxation import temperature_clock
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def reference_temperature(z, when, *, n0, ion_mass, length, speed):
    """Independent scalar integral along each exact advective characteristic."""
    unique, inverse = np.unique(z, return_inverse=True)
    te0 = 5.+.3*np.cos(2*np.pi*(unique-speed*when)/length)
    if when == 0:
        result = np.column_stack((te0, 6.-te0))
    else:
        result = np.array([temperature_clock(n0, ion_mass, t, 6.-t)(when) for t in te0])
    return result[inverse].reshape(np.shape(z)+(2,))


def run(destination, *, points=16, steps=80):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Fresh output directory required')
    if points not in (16,32) or type(steps) is not int or steps < 1:
        raise ValueError('Explicit 16/32 axial points and positive step count required')
    started = time.monotonic()
    identity = provenance(['examples.advected_coulomb_contact'])
    destination.mkdir(parents=True); freeze(destination, identity)
    n0=1e20; mi=39.948*atomic_mass; radius=.035; length=.24; speed=2400.; bias=.005
    duration=length/speed; dt=duration/steps
    fluid=StrongWallFluid(radius,length,4,points,mi+ME,degree=5,gammas=(5/3,5/3))
    mesh=MHDMesh(radius,length/2,8,points)
    original_geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    original=RegularPlasma(RegularFieldMap(fluid,original_geometry),mi)
    geometry=StrongWallGeometry(original_geometry)
    system=CoulombExchangeWallPlasma(RegularFieldMap(fluid,geometry),mi,
        resistivity=0.,electron_heat_fraction=1.)
    system.original_system=original
    state=np.zeros(fluid.shape);state[...,0]=n0;state[...,3]=n0*speed
    te0=5.+.3*np.cos(2*np.pi*fluid.z/length)
    state[...,4:]=n0*np.log(np.stack((te0,6.-te0),axis=-1))
    initial_state=state.copy()
    field=HarmonicMagnetic(mesh,np.zeros((mesh.nr+1,mesh.nz+1)),np.zeros((mesh.nr,mesh.nz)),bias)
    solver=GuardedModalGauss(system,state,field,absolute_tolerance=1e-12)
    sample=TransportSampler(fluid,mesh,20)
    volume=float(np.sum(sample.weight));particles=n0*volume
    temperature_initial_reference=reference_temperature(sample.Z,0.,n0=n0,ion_mass=mi,length=length,speed=speed)
    record=dict(terminal=False,success=False,steps=steps,dt=dt,duration=duration,
        completed_steps=0,spec=dict(elements=4,points=points,field_nr=8,field_nz=points),
        density_m3=n0,ion_mass_kg=mi,charge_state=1,gammas=[5/3,5/3],
        radius_m=radius,length_m=length,axial_speed_m_s=speed,bias_t=bias,
        initial_temperature='Te=5+0.3*cos(2*pi*z/L), Ti=6-Te [eV]',
        resistivity_ohm_m=0.,source_sha256=identity['source_sha256'],records=[],
        scope='Actual full moving/field equations with thermal exchange, no conduction or ionization. Independent exact reduced-model solution: relaxation along z-vt, constant n/v/total pressure and axial B. This is a component benchmark, not a complete collisional laboratory model or a cure for the separate magnetic/current accuracy problem.')

    def observe(u,b,when):
        n,v,_=sample.fluid(u)
        entropy=np.stack([sample.scalar(u[...,k])/n for k in (4,5)],axis=-1)
        temperature=np.exp(entropy+(2/3)*np.log(n[...,None]/n0))
        target=reference_temperature(sample.Z,when,n0=n0,ion_mass=mi,length=length,speed=speed)
        model=RegularCurrent(system.mapping,u[...,0],mi)
        packed=b.packed_increment(geometry);current=model.current(packed)
        j=sample.electric(current);increment=sample.magnetic(packed)
        magnetic=increment.copy();magnetic[...,2]+=bias
        thermal=1.5*E*n*np.sum(temperature,axis=-1)
        energy=float(np.sum(sample.weight*(thermal+.5*fluid.mass*n*np.sum(v*v,axis=-1)
            +np.sum(magnetic*magnetic,axis=-1)/(2*MU0)
            +model.reduced_mass*np.sum(j*j,axis=-1)/(2*n*E**2))))
        native=float(np.sum(fluid.volume*fluid.fields(u)['energy'])+.5*packed@(model.hessian@packed)
            +np.sum(mesh.volumes*model.current_energy(current))+bias*bias*float(np.sum(mesh.volumes))/(2*MU0))
        norm=lambda a:float(np.sqrt(np.sum(sample.weight*a*a)/volume))
        error=temperature-target;vd=v.copy();vd[...,2]-=speed
        pressure=n*E*np.sum(temperature,axis=-1);p0=n0*E*6.
        entropy_change=float(np.sum(fluid.volume[...,None]*(u[...,4:]-initial_state[...,4:]))/(particles*(2/3)))
        reference_entropy=float(1.5*np.sum(sample.weight*np.sum(np.log(target/temperature_initial_reference),axis=-1))/volume)
        return dict(time=when,temperature_rms_error_ev=[norm(error[...,k]) for k in (0,1)],
            temperature_max_error_ev=np.max(abs(error),axis=(0,1)).tolist(),
            temperature_min_ev=temperature.min(axis=(0,1)).tolist(),temperature_max_ev=temperature.max(axis=(0,1)).tolist(),
            density_relative_rms_error=norm((n-n0)/n0),pressure_relative_rms_error=norm((pressure-p0)/p0),
            maximum_velocity_error_m_s=float(np.linalg.norm(vd,axis=-1).max()),
            maximum_current_a_m2=float(np.linalg.norm(j,axis=-1).max()),
            maximum_field_departure_t=float(np.linalg.norm(increment,axis=-1).max()),
            maximum_drift_over_thermal_speed=float(np.max(np.linalg.norm(j,axis=-1)/(n*E*np.sqrt(E*(temperature[...,0]/ME+temperature[...,1]/mi))))),
            physical_energy_j=energy,native_energy_j=native,
            particle_relative_change=float(np.sum(sample.weight*sample.scalar(u[...,0]-initial_state[...,0]))/particles),
            entropy_change_per_initial_particle=entropy_change,
            reference_entropy_change_per_initial_particle=reference_entropy,
            entropy_change_error_per_initial_particle=entropy_change-reference_entropy)

    def save():
        record.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.monotonic()-started)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')

    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 4*1024**3 or shutil.disk_usage('.').free < 256*1024**2:
            raise RuntimeError('Moving relaxation exceeded own 4GiB/256MiB resource limits')

    np.savez_compressed(destination/'initial.npz',state=state,magnetic=field.packed_increment(geometry))
    record['initial_sha256']=sha(destination/'initial.npz');initial=observe(state,field,0.)
    record['initial']=initial;resources();save()
    for step in range(1,steps+1):
        try:
            updated,evolved,diagnostic=solver.step(state,field,dt)
            resources()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error))
            break
        raw=destination/f'step_{step:04}.npz'
        np.savez_compressed(raw,state=updated,magnetic=evolved.packed_increment(geometry),stage_changes=diagnostic.pop('stage_changes'))
        state,field=updated,evolved
        measured=observe(state,field,step*dt)
        measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,
            physical_energy_relative_change=(measured['physical_energy_j']-initial['physical_energy_j'])/initial['physical_energy_j'],
            native_energy_relative_change=(measured['native_energy_j']-initial['native_energy_j'])/initial['native_energy_j'])
        record['records'].append(measured);record['completed_steps']=step;resources();save()
        if step==1 or step%10==0:
            print('advected exchange',points,step,'/',steps,'T error',measured['temperature_rms_error_ev'],flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps,final_time=record['completed_steps']*dt)
    if provenance(['examples.advected_coulomb_contact']) != identity:
        raise RuntimeError('Source changed during moving relaxation')
    save();print('advected exchange terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--points',type=int,default=16)
    parser.add_argument('--steps',type=int,default=80);args=parser.parse_args();run(args.destination,points=args.points,steps=args.steps)
