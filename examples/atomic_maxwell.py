"""Complete vacuum and finite-mass electromagnetic wave evolutions."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, c as LIGHT_SPEED, elementary_charge as EV, epsilon_0 as EPS0
from scipy.optimize import brentq
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_maxwell import AtomicMaxwell
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance, freeze


def wave_reference(model, kind, amplitude, time):
    """Continuum linear mode; no discrete curl, force or timestep enters it."""
    k = np.pi/model.length
    mass = model.source.mass[[1, -1]]
    charge = EV*np.array([1., -1.])
    plasma2 = model.density*charge**2/(EPS0*mass) if kind != 'vacuum' else np.zeros(2)
    gyro = charge*model.axial_magnetic/mass
    if kind == 'circular':
        def dispersion(w):
            return w*w-(LIGHT_SPEED*k)**2-np.sum(plasma2*w/(w-gyro))
        lower = max(LIGHT_SPEED*k, gyro.max())*(1+1e-12)
        upper = 2*(LIGHT_SPEED*k+np.sqrt(plasma2.sum())+abs(gyro).max())
        frequency = brentq(dispersion, lower, upper, xtol=1e-4, rtol=1e-14)
    else:
        frequency = np.sqrt((LIGHT_SPEED*k)**2+plasma2.sum())
    phase = frequency*time
    face_sine = np.sin(k*np.arange(model.cells+1)*model.dx)
    face_sine[[0, -1]] = 0.  # Exact analytic conducting-end values.
    sine = np.sin(k*model.x)*np.sinc(.5/model.cells)
    cosine = np.cos(k*model.x)*np.sinc(.5/model.cells)
    electric = np.zeros((model.cells+1, 2))
    magnetic = np.zeros((model.cells, 2))
    velocity = np.zeros((model.cells, 3, 3))
    electric[:, 0] = amplitude*face_sine*np.cos(phase)
    if kind == 'circular':
        electric[:, 1] = -amplitude*face_sine*np.sin(phase)
        magnetic[:, 0] = k*amplitude/frequency*cosine*np.cos(phase)
        magnetic[:, 1] = -k*amplitude/frequency*cosine*np.sin(phase)
        speed = charge*amplitude/(mass*(frequency-gyro))
        velocity[:, [1, 2], 0] = sine[:, None]*speed*np.sin(phase)
        velocity[:, [1, 2], 1] = sine[:, None]*speed*np.cos(phase)
    else:
        magnetic[:, 1] = -k*amplitude/frequency*cosine*np.sin(phase)
        if kind != 'vacuum':
            velocity[:, [1, 2], 0] = sine[:, None]*(charge*amplitude/(mass*frequency))*np.sin(phase)
    current = model.density*np.einsum('s,csi->ci', charge, velocity[:, [1, 2]])
    return frequency, electric, magnetic, velocity, current


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh Maxwell evolution output required')
    rr = RadiativeRecombination(dict(nuclear_charge=1, transitions=[dict(initial_charge=1, kind='power_law',
        A_cm3_s=1e-12, reference_temperature_ev=2., exponent=0., minimum_temperature_ev=None, maximum_temperature_ev=None)]))
    network = AtomicNetwork(ConstantIonization(), rr, ionization_enabled=False, recombination_enabled=False, radiation=False)
    source = FixedFrameAtomicMoments(network, 39.948*atomic_mass, [200., -100., 300.])
    identity = provenance(['examples.atomic_maxwell'])
    out.mkdir(parents=True); freeze(out, identity)
    began = time.monotonic()
    report = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[], reference_density_m3=1e16,
        length_m=.1, temperature_ev=[.2, 2.], reference_velocity_m_s=source.reference_velocity.tolist(),
        scope='Full nonlinear spatial species/Maxwell waves, physical masses and light speed. Closed conducting verification interval; no reaction, free-boundary or formation claim.')
    histories = {}

    def save():
        report['elapsed_seconds'] = time.monotonic()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')

    save()
    for kind in ('vacuum', 'linear', 'circular'):
        for amplitude in ((1000.,) if kind == 'vacuum' else (100., 1000.)):
            for cells in (16, 32, 64):
                model = AtomicMaxwell(source, cells, .1, 1e16, axial_magnetic=.02 if kind == 'circular' else 0.,
                                      maximum_plasma_phase=.04*16/cells)
                frequency, electric, magnetic, velocity, _ = wave_reference(model, kind, amplitude, 0.)
                atom = network.initial([0., 0. if kind == 'vacuum' else 1.], 2.)
                state = np.asarray([source.initial(atom, .2, v) for v in velocity])
                duration = 2.3*2*np.pi/frequency
                times = [0.]; states = [state.copy()]; fields = [electric.copy()]; magnets = [magnetic.copy()]
                t = 0.; failure = None; max_work = 0.
                while t < duration:
                    dt = min(model.timestep(state, electric, magnetic), duration-t)
                    try:
                        candidate = model.step(state, electric, magnetic, t, dt, reactions=False)
                    except (ValueError, RuntimeError) as error:
                        failure = dict(step=len(times), time_s=t, dt_s=dt, kind=type(error).__name__, message=str(error))
                        break
                    state, electric, magnetic = candidate
                    t = min(t+dt, duration)
                    times.append(t); states.append(state.copy()); fields.append(electric.copy()); magnets.append(magnetic.copy())
                    max_work = max(max_work, abs(model.last_diagnostic['transverse_work_residual_ev_m_s']))
                times = np.asarray(times); states = np.asarray(states); fields = np.asarray(fields); magnets = np.asarray(magnets)
                references = [wave_reference(model, kind, amplitude, t) for t in times]
                exact_e = np.asarray([r[1] for r in references]); exact_b = np.asarray([r[2] for r in references])
                exact_j = np.asarray([r[4] for r in references])
                currents = np.asarray([EV*model.density*np.einsum('s,csi->ci', model.charge, model.number_velocity(s)) for s in states])
                rows = [model.totals(s, e, b) for s, e, b in zip(states, fields, magnets)]
                budgets = {key: np.asarray([r[key] for r in rows]) for key in rows[0]}
                electric_z = np.asarray([model.electric(s) for s in states])
                initial_wave_energy = budgets['transverse_electric_energy_ev'][0]+budgets['transverse_magnetic_energy_ev'][0]
                error_j = np.linalg.norm(currents-exact_j, axis=2).max()
                case = dict(name=f'{kind}_e{amplitude:g}_n{cells}', kind=kind, cells=cells, electric_amplitude_v_m=amplitude,
                    axial_magnetic_t=model.axial_magnetic, frequency_rad_s=float(frequency), duration_s=duration,
                    steps=len(times)-1, success=failure is None, failure=failure, completed_time_s=float(times[-1]),
                    maximum_dt_s=float(np.diff(times).max(initial=0.)),
                    maximum_electric_error_over_amplitude=float(np.linalg.norm(fields-exact_e, axis=2).max()/amplitude),
                    maximum_magnetic_error_over_amplitude=float(np.linalg.norm(magnets-exact_b, axis=2).max()/(np.pi/model.length*amplitude/frequency)),
                    maximum_current_error_a_m2=float(error_j),
                    maximum_current_error_over_mode_amplitude=float(error_j/np.linalg.norm(exact_j, axis=2).max()) if kind != 'vacuum' else None,
                    maximum_longitudinal_electric_v_m=float(abs(electric_z).max()),
                    maximum_nucleus_integral_error=float(abs(budgets['nuclei']-budgets['nuclei'][0]).max()),
                    maximum_charge_integral_error=float(abs(budgets['charge']-budgets['charge'][0]).max()),
                    maximum_total_energy_error_ev_m=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()),
                    maximum_energy_error_over_initial_transverse_field=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()/initial_wave_energy),
                    maximum_linear_field_energy_error_over_initial_transverse_field=float(abs(budgets['linear_field_total_energy_ev']-budgets['linear_field_total_energy_ev'][0]).max()/initial_wave_energy),
                    maximum_gauss_residual_over_reference_charge=float(budgets['maximum_gauss_residual_over_reference_charge'].max()),
                    maximum_stage_transverse_work_residual_ev_m_s=float(max_work),
                    maximum_neutral_number=float(states[:, :, 0].max()),
                    maximum_total_species_number=float(states[:, :, model.number].sum(axis=2).max()))
                history = dict(time_s=times, states=states, transverse_electric_v_m=fields, transverse_magnetic_t=magnets,
                    longitudinal_electric_v_m=electric_z, current_a_m2=currents, analytic_electric_v_m=exact_e,
                    analytic_magnetic_t=exact_b, analytic_current_a_m2=exact_j, x_m=model.x, face_x_m=np.arange(cells+1)*model.dx, **budgets)
                np.savez_compressed(out/(case['name']+'.npz'), **history)
                histories[kind, amplitude, cells] = history
                report['cases'].append(case); save(); print('Spatial Maxwell wave', case, flush=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for cells in (16, 32, 64):
        h = histories['circular', 100., cells]
        axes[0, 0].plot(h['face_x_m'], h['transverse_electric_v_m'][-1, :, 1], label=f'{cells} cells')
        axes[0, 1].plot(h['x_m'], h['transverse_magnetic_t'][-1, :, 1]*1e6, label=f'{cells} cells')
    axes[0, 0].plot(h['face_x_m'], h['analytic_electric_v_m'][-1, :, 1], 'k--', label='Continuum circular mode')
    axes[0, 1].plot(h['x_m'], h['analytic_magnetic_t'][-1, :, 1]*1e6, 'k--', label='Continuum circular mode')
    for kind in ('vacuum', 'linear', 'circular'):
        amplitude = 1000. if kind == 'vacuum' else 100.
        rows = [r for r in report['cases'] if r['kind'] == kind and r['electric_amplitude_v_m'] == amplitude]
        axes[1, 0].loglog([r['cells'] for r in rows], [r['maximum_electric_error_over_amplitude'] for r in rows], 'o-', label=kind)
        h = histories[kind, amplitude, 64]
        scale = h['transverse_electric_energy_ev'][0]+h['transverse_magnetic_energy_ev'][0]
        axes[1, 1].plot(h['time_s']*1e9, (h['total_energy_ev']-h['total_energy_ev'][0])/scale, label=kind)
    for ax, xlabel, ylabel in zip(axes.flat, ('Axial coordinate [m]', 'Axial coordinate [m]', 'Cells', 'Time [ns]'),
            ('Final Ey [V/m]', 'Final By [µT]', 'Maximum vector E error / amplitude', 'Energy drift / initial transverse field energy')):
        ax.set(xlabel=xlabel, ylabel=ylabel); ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Spatial Maxwell / physical-mass species waves\nConducting verification interval; 20 mT bias in circular case')
    fig.savefig(out/'spatial_maxwell.png', dpi=170); fig.savefig(out/'spatial_maxwell.pdf'); plt.close(fig)
    if provenance(['examples.atomic_maxwell']) != identity:
        raise RuntimeError('Source changed during spatial Maxwell evolution')
    report['terminal'] = True; report['success'] = all(r['success'] for r in report['cases']); save()
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('destination'); a = p.parse_args()
    r = run(a.destination); raise SystemExit(0 if r['success'] else 1)
