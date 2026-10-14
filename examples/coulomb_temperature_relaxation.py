"""Complete moving-plasma relaxation against an independent scalar clock."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0, mu_0 as MU0
from scipy.integrate import quad
from scipy.optimize import brentq
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
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def temperature_clock(n, mi, initial_te, initial_ti):
    """Independent separated scalar integral; no production collision call."""
    total = initial_te+initial_ti
    mean = total/2
    coefficient = 8*np.sqrt(2*np.pi)*n*(E**2/(4*np.pi*EPS0))**2/(3*ME*mi*E**1.5)

    upper = np.log(initial_te-mean)

    def elapsed(log_difference):
        def integrand(y):
            t = mean+np.exp(y)
            log = 23.-.5*np.log(n/1e6)+1.5*np.log(t)
            return (t/ME+(total-t)/mi)**1.5/(2*coefficient*log)
        return quad(integrand, log_difference, upper, epsabs=2e-14, epsrel=2e-11)[0]

    def expected(when):
        if when == 0:
            return np.array([initial_te, initial_ti])
        y = brentq(lambda y: elapsed(y)-when, np.log(1e-8), upper, xtol=1e-13)
        te = mean+np.exp(y)
        return np.array([te, total-te])

    return expected


def run(destination, *, steps=40):
    destination = Path(destination)
    if destination.exists() or type(steps) is not int or steps < 1:
        raise ValueError('Fresh output path and positive integer steps required')
    started = time.monotonic()
    identity = provenance(['examples.coulomb_temperature_relaxation'])
    destination.mkdir(parents=True); freeze(destination, identity)
    n0 = 1e20; mi = 39.948*atomic_mass; initial_temperature = np.array([5., 1.])
    duration = 4e-4; dt = duration/steps
    fluid = StrongWallFluid(.035, .24, 4, 16, mi+ME, degree=5, gammas=(5/3, 5/3))
    mesh = MHDMesh(.035, .12, 8, 16)
    original_geometry = CompatibleGeometry(mesh, degree=5, quadrature_order=12)
    original = RegularPlasma(RegularFieldMap(fluid, original_geometry), mi)
    geometry = StrongWallGeometry(original_geometry)
    system = CoulombExchangeWallPlasma(RegularFieldMap(fluid, geometry), mi,
        resistivity=0., electron_heat_fraction=1.)
    system.original_system = original
    state = np.zeros(fluid.shape); state[..., 0] = n0
    state[..., 4:] = n0*np.log(initial_temperature)
    initial_state = state.copy()
    field = HarmonicMagnetic(mesh, np.zeros((mesh.nr+1, mesh.nz+1)), np.zeros((mesh.nr, mesh.nz)), 0.)
    solver = GuardedModalGauss(system, state, field, absolute_tolerance=1e-12)
    sample = TransportSampler(fluid, mesh, 20)
    expected = temperature_clock(n0, mi, *initial_temperature)
    initial_energy = float(np.sum(fluid.volume*fluid.fields(state)['energy']))
    particles = float(np.sum(fluid.volume*n0))
    record = dict(terminal=False, success=False, steps=steps, dt=dt, duration=duration,
        completed_steps=0, density_m3=n0, ion_mass_kg=mi, charge_state=1,
        initial_temperature_ev=initial_temperature.tolist(), equilibrium_temperature_ev=3.,
        gammas=[5/3,5/3], resistivity_ohm_m=0., bias_t=0.,
        spec=dict(elements=4, points=16, field_nr=8, field_nz=16),
        source_sha256=identity['source_sha256'], records=[],
        scope='Full moving/field equations with conditional Maxwellian thermal exchange. Uniform zero-current, zero-drift fixture; scalar quadrature baseline. No numerical damping, vacuum, formation, neutral or ionization model.')
    # Independent rounded NRL cgs/eV coefficient: masses in grams, n in cm^-3.
    log0 = 23.-np.log(np.sqrt(n0/1e6)*initial_temperature[0]**-1.5)
    nrl = 1.8e-19*np.sqrt((ME*1000)*(mi*1000))*(n0/1e6)*log0
    nrl /= (ME*1000*initial_temperature[1]+mi*1000*initial_temperature[0])**1.5
    first_rate = system.derivative(state, field)
    actual_frequency = first_rate.diagnostics['exchange_frequency_range_s'][0]
    record.update(initial_exchange_rate_s=actual_frequency, nrl_rounded_initial_rate_s=float(nrl),
        rate_relative_difference_from_rounded_nrl=float(actual_frequency/nrl-1),
        initial_exchange_power_to_electrons=first_rate.diagnostics['exchange_power_to_electrons'])

    def save():
        record.update(elapsed_seconds=time.monotonic()-started,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')

    def resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 1024**3 or shutil.disk_usage('.').free < 256*1024**2:
            raise RuntimeError('Relaxation example exceeded own 1GiB/256MiB resource limits')

    np.savez_compressed(destination/'initial.npz', state=state, magnetic=field.packed_increment(geometry))
    record['initial_sha256'] = sha(destination/'initial.npz'); resources(); save()
    for step in range(1, steps+1):
        try:
            updated, evolved, diagnostic = solver.step(state, field, dt)
            resources()
        except Exception as error:
            record.update(failed_step=step, failure_type=type(error).__name__, failure_message=str(error))
            break
        raw = destination/f'step_{step:04}.npz'
        packed = evolved.packed_increment(geometry)
        np.savez_compressed(raw, state=updated, magnetic=packed, stage_changes=diagnostic.pop('stage_changes'))
        state, field = updated, evolved
        n, velocity, _ = sample.fluid(state)
        entropy = np.stack([sample.scalar(state[...,k])/n for k in (4,5)], axis=-1)
        temperature = np.exp(entropy+(2/3)*np.log(n[...,None]/n0))
        target = expected(step*dt)
        model = RegularCurrent(system.mapping, state[...,0], mi)
        current = model.current(packed)
        b = sample.magnetic(packed); j = sample.electric(current)
        thermal = 1.5*E*n*np.sum(temperature, axis=-1)
        # No harmonic bias in this uniform relaxation fixture.
        physical = float(np.sum(sample.weight*(thermal+.5*fluid.mass*n*np.sum(velocity**2, axis=-1)
            +np.sum(b*b, axis=-1)/(2*MU0)+model.reduced_mass*np.sum(j*j, axis=-1)/(2*n*E**2))))
        native = float(np.sum(fluid.volume*fluid.fields(state)['energy'])+.5*packed@(model.hessian@packed)
            +np.sum(mesh.volumes*model.current_energy(current)))
        record['records'].append(dict(step=step,time=step*dt,arrays_file=raw.name,arrays_sha256=sha(raw),solver=diagnostic,
            expected_temperature_ev=target.tolist(),
            temperature_mean_ev=(np.sum(sample.weight[...,None]*temperature,axis=(0,1))/np.sum(sample.weight)).tolist(),
            maximum_temperature_error_ev=np.max(abs(temperature-target),axis=(0,1)).tolist(),
            physical_energy_relative_drift=(physical-initial_energy)/initial_energy,
            native_energy_relative_drift=(native-initial_energy)/initial_energy,
            particle_relative_drift=float(np.sum(sample.weight*(n-n0))/particles),
            entropy_change_per_initial_particle=float(np.sum(fluid.volume[...,None]*(state[...,4:]-initial_state[...,4:]))/(particles*(2/3))),
            minimum_density=float(n.min()), maximum_speed=float(np.linalg.norm(velocity,axis=-1).max()),
            maximum_current=float(np.linalg.norm(j,axis=-1).max()),maximum_field=float(np.linalg.norm(b,axis=-1).max())))
        record['completed_steps'] = step; save()
        if step == 1 or step%10 == 0:
            print('Coulomb relaxation',step,'/',steps,'T',record['records'][-1]['temperature_mean_ev'],flush=True)
    record.update(terminal=True, success=record['completed_steps']==steps)
    if provenance(['examples.coulomb_temperature_relaxation']) != identity:
        raise RuntimeError('Source changed during relaxation')
    save(); print('Coulomb relaxation terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--steps',type=int,default=40)
    args=parser.parse_args();run(args.destination,steps=args.steps)
