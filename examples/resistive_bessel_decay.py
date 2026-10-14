"""Full moving, heated plasma compared with independent linear Bessel decay."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.special import j0, j1, jn_zeros
from plasma_bench.guarded_modal_gauss import GuardedModalGauss
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.resistive_wall_plasma import ResistiveWallPlasma
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha
from .strong_wall_reference import initialize


def run(destination, *, nr=8, steps=40, amplitude=1e-4, resistivity=1e-4):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Fresh output directory required')
    if nr not in (8, 16) or type(steps) is not int or steps < 1:
        raise ValueError('This small experiment supports nr8/16 and positive steps')
    if not np.isfinite([amplitude, resistivity]).all() or min(amplitude, resistivity) <= 0:
        raise ValueError('Positive finite magnetic amplitude and resistivity required')
    started = time.monotonic()
    identity = provenance(['examples.resistive_bessel_decay'])
    destination.mkdir(parents=True)
    freeze(destination, identity)
    spec = dict(kind='contact', elements=nr//2, points=16, field_nr=nr, field_nz=16)
    ref, fluid, _, original, _ = initialize(spec)
    system = ResistiveWallPlasma(original.mapping, ref.ion_mass, resistivity=resistivity,
        electron_heat_fraction=1., **original.solver_options)
    system.original_system = original.original_system
    n0 = ref.number_density
    state = np.zeros(fluid.shape)
    state[..., 0] = n0
    state[..., 4:] = n0*np.log(np.array([ref.electron_ev, ref.ion_ev]))
    mesh = system.geometry.mesh
    q = float(jn_zeros(1, 1)[0]/ref.radius)
    r, z = mesh.corners
    psi = amplitude*r*j1(q*r)/q+np.zeros_like(z)
    psi[[0, -1]] = 0.
    field = HarmonicMagnetic(mesh, psi, np.zeros((mesh.nr, mesh.nz)), 0.)
    reduced_mass = ref.ion_mass*ref.electron_mass/(ref.ion_mass+ref.electron_mass)
    de2 = reduced_mass/(MU0*n0*E**2)
    decay_rate = resistivity*q*q/(MU0*(1+de2*q*q))
    duration = 5/decay_rate
    dt = duration/steps
    solver = GuardedModalGauss(system, state, field, absolute_tolerance=1e-12)
    sample = TransportSampler(fluid, mesh, 20)
    initial_bz = amplitude*j0(q*sample.R)+np.zeros_like(sample.Z)
    initial_jtheta = amplitude*q*j1(q*sample.R)/MU0+np.zeros_like(sample.Z)
    norm = lambda value: float(np.sqrt(np.sum(sample.weight*value*value)))
    signal = norm(initial_bz)
    analytic_energy = float(np.sum(sample.weight*(initial_bz**2/(2*MU0)
        +reduced_mass*initial_jtheta**2/(2*n0*E**2))))
    record = dict(terminal=False, success=False, spec=spec, amplitude_t=amplitude,
        resistivity_ohm_m=resistivity, electron_heat_fraction=1., density_m3=n0,
        temperatures_ev=[ref.electron_ev, ref.ion_ev], q_per_m=q,
        finite_mass_decay_rate_s=decay_rate, duration=duration, steps=steps, dt=dt,
        source_sha256=identity['source_sha256'], completed_steps=0, records=[],
        scope='Actual nonlinear moving resistive plasma; small-amplitude linear Bessel target, not exact nonlinear solution or formation. Electron-only prescribed heat partition; no temperature exchange or transport closure.')

    def model_for(u):
        return RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)

    def observe(u, b, when):
        model = model_for(u)
        packed = b.packed_increment(system.geometry)
        current = model.current(packed)
        magnetic = sample.magnetic(packed)
        j = sample.electric(current)
        n, velocity, _ = sample.fluid(u)
        entropy = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
        temperature = np.exp(entropy+(fluid.thermal.gammas-1)*np.log(n[..., None]/n0))
        thermal = n*E*np.sum(temperature/(fluid.thermal.gammas-1), axis=-1)
        kinetic = .5*fluid.mass*n*np.sum(velocity*velocity, axis=-1)
        magnetic_energy = np.sum(magnetic*magnetic, axis=-1)/(2*MU0)
        current_energy = reduced_mass*np.sum(j*j, axis=-1)/(2*n*E**2)
        physical = {name:float(np.sum(sample.weight*value)) for name,value in
            [('thermal', thermal), ('kinetic', kinetic), ('magnetic', magnetic_energy), ('current', current_energy)]}
        native = float(np.sum(fluid.volume*fluid.fields(u)['energy'])
            +.5*packed@(model.hessian@packed)+np.sum(mesh.volumes*model.current_energy(current)))
        expected = initial_bz*np.exp(-decay_rate*when)
        error = norm(magnetic[..., 2]-expected)
        return dict(time=when, physical_energy=physical, native_total_energy=native,
            bz_error_over_initial=error/signal,
            bz_error_over_remaining=error/(signal*np.exp(-decay_rate*when)),
            jtheta_error_over_initial=norm(j[..., 1]-initial_jtheta*np.exp(-decay_rate*when))/norm(initial_jtheta),
            ohmic_power=float(np.sum(sample.weight*resistivity*np.sum(j*j, axis=-1))),
            linear_released_energy=analytic_energy*(-np.expm1(-2*decay_rate*when)),
            minimum_density=float(n.min()), maximum_speed=float(np.sqrt(np.sum(velocity*velocity, axis=-1)).max()),
            minimum_temperature_ev=temperature.min(axis=(0,1)).tolist())

    def save():
        record.update(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.monotonic()-started)
        (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')

    def check_resources():
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 1024**3 or shutil.disk_usage('.').free < 256*1024**2:
            raise RuntimeError('Small experiment exceeded its 1GiB/256MiB resource budget')

    check_resources()
    np.savez_compressed(destination/'initial.npz', state=state, magnetic=field.packed_increment(system.geometry))
    initial = observe(state, field, 0.)
    record['initial'] = initial
    record['initial_sha256'] = sha(destination/'initial.npz')
    save()
    for step in range(1, steps+1):
        try:
            updated, evolved, diagnostic = solver.step(state, field, dt)
            check_resources()
        except Exception as error:
            record.update(failed_step=step, failure_type=type(error).__name__, failure_message=str(error))
            break
        raw = destination/f'step_{step:04}.npz'
        np.savez_compressed(raw, state=updated, magnetic=evolved.packed_increment(system.geometry),
            stage_changes=diagnostic.pop('stage_changes'))
        state, field = updated, evolved
        measured = observe(state, field, step*dt)
        measured.update(step=step, arrays_file=raw.name, arrays_sha256=sha(raw),
            solver=diagnostic,
            native_energy_drift_over_initial_field=(measured['native_total_energy']-initial['native_total_energy'])/analytic_energy,
            physical_energy_drift_over_initial_field=(sum(measured['physical_energy'].values())-sum(initial['physical_energy'].values()))/analytic_energy)
        record['records'].append(measured)
        record['completed_steps'] = step
        save()
        if step == 1 or step % 10 == 0:
            print('resistive decay', nr, step, '/', steps, 'B error', measured['bz_error_over_remaining'], flush=True)
    record.update(terminal=True, success=record['completed_steps']==steps)
    if provenance(['examples.resistive_bessel_decay']) != identity:
        raise RuntimeError('Source changed during resistive evolution')
    save()
    print('resistive decay terminal', record['success'], record.get('failure_message', ''), flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    parser.add_argument('--nr', type=int, default=8)
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--amplitude', type=float, default=1e-4)
    args = parser.parse_args()
    run(args.destination, nr=args.nr, steps=args.steps, amplitude=args.amplitude)
