"""Magnetic decay while collisional cooling changes electrical resistance."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0, mu_0 as MU0
from scipy.integrate import quad
from scipy.special import j0, j1, jn_zeros
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.guarded_modal_gauss import GuardedModalGauss
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_plasma import RegularPlasma
from plasma_bench.spitzer_wall_plasma import SpitzerWallPlasma
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.strong_wall_plasma import StrongWallFluid
from .coulomb_temperature_relaxation import temperature_clock
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def resistance_reference(te, n):
    """SI low-field Spitzer law, independent of the production coefficient code."""
    log = 23.-np.log(np.sqrt(n/1e6)/te**1.5)
    return .51*np.sqrt(ME)*E**2*log/(3*(2*np.pi)**1.5*EPS0**2*(E*te)**1.5)


def run(destination, *, nr=8, steps=80, amplitude=1e-4):
    destination = Path(destination)
    if destination.exists() or nr not in (8, 16) or type(steps) is not int or steps < 1:
        raise ValueError('Fresh output, nr8/16 and positive integer steps required')
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError('Positive finite field amplitude required')
    started = time.monotonic()
    identity = provenance(['examples.spitzer_bessel_decay'])
    destination.mkdir(parents=True); freeze(destination, identity)
    n0=1e21; mi=39.948*atomic_mass; radius=.035; length=.24
    duration=12e-6; dt=duration/steps
    fluid=StrongWallFluid(radius,length,nr//2,16,mi+ME,degree=5,gammas=(5/3,5/3))
    mesh=MHDMesh(radius,length/2,nr,16)
    original_geometry=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
    original=RegularPlasma(RegularFieldMap(fluid,original_geometry),mi)
    geometry=StrongWallGeometry(original_geometry)
    system=SpitzerWallPlasma(RegularFieldMap(fluid,geometry),mi)
    system.original_system=original
    state=np.zeros(fluid.shape); state[...,0]=n0
    state[...,4:]=n0*(np.log([5.,1.])-(fluid.thermal.gammas-1)*np.log(n0/fluid.thermal.n0))
    initial_state=state.copy()
    q=float(jn_zeros(1,1)[0]/radius); r,z=mesh.corners
    psi=amplitude*r*j1(q*r)/q+np.zeros_like(z); psi[[0,-1]]=0.
    field=HarmonicMagnetic(mesh,psi,np.zeros((mesh.nr,mesh.nz)),0.)
    reduced=ME*mi/(ME+mi); de2=reduced/(MU0*n0*E**2)
    solver=GuardedModalGauss(system,state,field,absolute_tolerance=1e-12)
    sample=TransportSampler(fluid,mesh,20); volume=float(np.sum(sample.weight))
    bz0=amplitude*j0(q*sample.R)+np.zeros_like(sample.Z)
    jtheta0=amplitude*q*j1(q*sample.R)/MU0+np.zeros_like(sample.Z)
    norm=lambda a:float(np.sqrt(np.sum(sample.weight*a*a)))
    signal=norm(bz0); jsignal=norm(jtheta0)
    signal_energy=float(np.sum(sample.weight*(bz0*bz0/(2*MU0)+reduced*jtheta0*jtheta0/(2*n0*E**2))))
    temperatures=temperature_clock(n0,mi,5.,1.)

    def linear_target(when):
        temperature=temperatures(when)
        integral=quad(lambda t:resistance_reference(temperatures(t)[0],n0),0.,when,
            epsabs=1e-22,epsrel=2e-10)[0] if when else 0.
        return temperature,float(np.exp(-q*q*integral/(MU0*(1+de2*q*q))))

    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,dt=dt,duration=duration,
        spec=dict(elements=nr//2,points=16,field_nr=nr,field_nz=16),radius_m=radius,length_m=length,
        density_m3=n0,ion_mass_kg=mi,gammas=[5/3,5/3],initial_temperature_ev=[5.,1.],
        amplitude_t=amplitude,q_per_m=q,initial_analytic_field_current_energy_j=signal_energy,
        source_sha256=identity['source_sha256'],records=[],
        scope='Full moving Z1 weak-field Spitzer friction/electron Joule heat plus thermal exchange. Independent time-dependent linear Bessel reference; nonlinear forces retained. No conduction/thermal forces/viscosity/neutrals/vacuum/formation.')

    def observe(u,b,when):
        model=RegularCurrent(system.mapping,u[...,0],mi); packed=b.packed_increment(geometry)
        current=model.current(packed); magnetic=sample.magnetic(packed); j=sample.electric(current)
        n,v,_=sample.fluid(u)
        entropy=np.stack([sample.scalar(u[...,k])/n for k in (4,5)],axis=-1)
        temperature=np.exp(entropy+(2/3)*np.log(n[...,None]/fluid.thermal.n0))
        expected,factor=linear_target(when)
        thermal=1.5*n*E*np.sum(temperature,axis=-1)
        values=dict(thermal=thermal,kinetic=.5*fluid.mass*n*np.sum(v*v,axis=-1),
            magnetic=np.sum(magnetic*magnetic,axis=-1)/(2*MU0),
            current=reduced*np.sum(j*j,axis=-1)/(2*n*E**2))
        energy={k:float(np.sum(sample.weight*a)) for k,a in values.items()}
        native=float(np.sum(fluid.volume*fluid.fields(u)['energy'])+.5*packed@(model.hessian@packed)
            +np.sum(mesh.volumes*model.current_energy(current)))
        diagnostic=system.derivative(u,b).diagnostics
        error=norm(magnetic[...,2]-bz0*factor)
        return dict(time=when,linear_decay_factor=factor,reference_temperature_ev=expected.tolist(),
            reference_resistivity_ohm_m=float(resistance_reference(expected[0],n0)),
            bz_error_over_initial=error/signal,bz_error_over_remaining=error/(signal*factor),
            jtheta_error_over_initial=norm(j[...,1]-jtheta0*factor)/jsignal,
            jtheta_error_over_remaining=norm(j[...,1]-jtheta0*factor)/(jsignal*factor),
            bz_mode_amplitude=float(np.sum(sample.weight*magnetic[...,2]*bz0)/signal**2),
            temperature_mean_ev=(np.sum(sample.weight[...,None]*temperature,axis=(0,1))/volume).tolist(),
            temperature_rms_departure_ev=[norm(temperature[...,k]-expected[k])/np.sqrt(volume) for k in (0,1)],
            physical_energy_j=energy,native_total_energy_j=native,
            linear_released_energy_j=signal_energy*(1-factor*factor),
            maximum_speed_m_s=float(np.linalg.norm(v,axis=-1).max()),
            minimum_density_m3=float(n.min()),minimum_temperature_ev=temperature.min(axis=(0,1)).tolist(),
            particle_relative_change=float(np.sum(sample.weight*sample.scalar(u[...,0]-initial_state[...,0]))/(n0*volume)),
            transport={k:diagnostic[k] for k in ('ohmic_power','resistivity_range_ohm_m',
                'electron_collision_rate_range_s','spitzer_coulomb_log_range','maximum_electron_hall_parameter',
                'maximum_drift_over_species_thermal_speed','weak_field_review','relative_drift_review')},
            instantaneous_native_power_sum_w=diagnostic['total_power'])

    def save():
        record.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.monotonic()-started)
        (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')

    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 4*1024**3 or shutil.disk_usage('.').free < 256*1024**2:
            raise RuntimeError('Spitzer decay exceeded own4GiB/256MiB limits')

    np.savez_compressed(destination/'initial.npz',state=state,magnetic=field.packed_increment(geometry))
    record['initial_sha256']=sha(destination/'initial.npz'); initial=observe(state,field,0.)
    record['initial']=initial; resources(); save()
    for step in range(1,steps+1):
        try:
            updated,evolved,diagnostic=solver.step(state,field,dt); resources()
        except Exception as error:
            record.update(failed_step=step,failure_type=type(error).__name__,failure_message=str(error)); break
        raw=destination/f'step_{step:04}.npz'
        np.savez_compressed(raw,state=updated,magnetic=evolved.packed_increment(geometry),stage_changes=diagnostic.pop('stage_changes'))
        state,field=updated,evolved; measured=observe(state,field,step*dt)
        measured.update(step=step,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,
            physical_energy_drift_over_initial_field=(sum(measured['physical_energy_j'].values())-sum(initial['physical_energy_j'].values()))/signal_energy,
            native_energy_drift_over_initial_field=(measured['native_total_energy_j']-initial['native_total_energy_j'])/signal_energy)
        record['records'].append(measured); record['completed_steps']=step; save()
        if step==1 or step%10==0:
            print('Spitzer decay',nr,step,'/',steps,'B error',measured['bz_error_over_remaining'],flush=True)
    record.update(terminal=True,success=record['completed_steps']==steps)
    if provenance(['examples.spitzer_bessel_decay']) != identity:
        raise RuntimeError('Source changed during Spitzer decay')
    save(); print('Spitzer decay terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--nr',type=int,default=8); parser.add_argument('--steps',type=int,default=80)
    parser.add_argument('--amplitude',type=float,default=1e-4)
    args=parser.parse_args(); run(args.destination,nr=args.nr,steps=args.steps,amplitude=args.amplitude)
