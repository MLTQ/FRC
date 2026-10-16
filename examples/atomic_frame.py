"""Full heated argon histories in fixed velocity coordinates."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from .atomic_moments import evolve
from .atomic_network import argon_models
from .implicit_plasma_records import provenance, freeze


def run(destination, baseline, neutral_path, ion_path, ionization_path):
    out = Path(destination); baseline = Path(baseline)
    if out.exists():raise ValueError('Fresh output directory required')
    baseline_report = json.loads((baseline/'summary.json').read_text())
    reference_case = next(c for c in baseline_report['cases'] if c['name']=='heated_seed_1400')
    if not baseline_report['terminal'] or not reference_case['success']:
        raise ValueError('A terminal complete unboosted reference history is required')
    with np.load(baseline/'heated_seed_1400.npz') as saved:
        reference = {k:saved[k].copy() for k in saved.files}
    models, inputs = argon_models(neutral_path, ion_path, ionization_path)
    network = models['joint']; density = 8.961429541213399e19
    velocity = np.zeros((network.count+1, 3))
    velocity[:-1, 2] = np.where(network.charge==0, 1000., -1000.)
    velocity[-1] = [0., 10000., 0.]
    boost = np.array([20000., -10000., 5000.])
    atomic = network.initial(np.r_[.9, .1, np.zeros(17)], 10.)
    identity = provenance(['examples.atomic_frame']); out.mkdir(parents=True); freeze(out, identity)
    report = dict(terminal=False, source_sha256=identity['source_sha256'], inputs=inputs,
        baseline=str(baseline.resolve()), baseline_source_sha256=baseline_report['source_sha256'],
        density_m3=density, common_physical_boost_m_s=boost.tolist(), cases=[],
        scope='Full 3.5 us heated argon reaction source; fixed Galilean coordinates. No electromagnetic or spatial feedback; existing small-drift and photon/atomic-data limits remain.')
    started = time.monotonic(); histories = {}
    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    for name, frame in [('boost_frame', boost),
                        ('offset_frame', boost+np.array([137., -211., 93.]))]:
        model = FixedFrameAtomicMoments(network, 39.948*atomic_mass, frame)
        initial = model.initial(atomic, .1, velocity+boost)
        summary, history = evolve(model, initial, density, 3.5e-6, 1400, 20.)
        history['laboratory_states'] = np.asarray([model.laboratory_state(s) for s in history['states']])
        history['reference_velocity_m_s'] = frame
        summary.update(name=name, reference_velocity_m_s=frame.tolist())
        if summary['success']:
            expected_momentum = reference['states'][:, model.momentum].copy()
            numbers = np.column_stack((reference['states'][:, :network.count], reference['electrons']))
            expected_momentum += (model.mass[:, None]/model.momentum_scale*boost)*numbers[:, :, None]
            expected_energy = reference['states'][:, model.energy].copy()
            expected_energy += (model.momentum_scale/EV)*(reference['states'][:, model.momentum]@boost)
            expected_energy += (.5*(model.mass/EV)*(boost@boost))*numbers
            expected_current = reference['current_a_m2']-EV*density*reference['charge_balance'][:, None]*boost
            summary['reference_comparison'] = dict(
                maximum_population_l1=float(abs(history['states'][:, :network.count]-reference['states'][:, :network.count]).sum(axis=1).max()),
                maximum_temperature_difference_ev=float(abs(history['electron_temperature_ev']-reference['electron_temperature_ev']).max()),
                maximum_current_difference_a_m2=float(np.linalg.norm(history['current_a_m2']-expected_current, axis=1).max()),
                maximum_species_lab_momentum_difference_scaled=float(abs(history['laboratory_states'][:, model.momentum]-expected_momentum).max()),
                maximum_species_lab_energy_difference_ev=float(abs(history['laboratory_states'][:, model.energy]-expected_energy).max()))
        np.savez_compressed(out/f'{name}.npz', **history)
        histories[name] = history; report['cases'].append(summary); save()
        print('Atomic frame', name, summary['success'], summary.get('reference_comparison'), flush=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for name, history in histories.items():
        t = history['time_s']*1e6; end = len(t)
        axes[0, 0].plot(t, history['electrons'], label=name)
        axes[0, 1].plot(t, history['electron_temperature_ev'], label=name)
        axes[1, 0].plot(t, history['electron_temperature_ev']-reference['electron_temperature_ev'][:end], label=name)
        axes[1, 1].plot(t, history['energy_budget_ev']-history['energy_budget_ev'][0], label=name)
    for ax, label in zip(axes.flat, ('Electrons / nucleus', 'Electron temperature [eV]',
            'Temperature minus unboosted reference [eV]', 'Lab energy budget residual [eV/nucleus]')):
        ax.set(xlabel='Time [µs]', ylabel=label); ax.grid(alpha=.2); ax.legend()
    fig.suptitle('Same heated argon source with common velocity boost; fixed-frame coordinates')
    fig.savefig(out/'atomic_frames.png', dpi=170); fig.savefig(out/'atomic_frames.pdf'); plt.close(fig)
    if provenance(['examples.atomic_frame'])!=identity:raise RuntimeError('Source changed during evolution')
    report['terminal'] = True; report['success'] = all(c['success'] for c in report['cases']); save()
    print('Atomic frame terminal', report['success'], out, flush=True)
    return report


if __name__=='__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('baseline')
    parser.add_argument('neutral_input'); parser.add_argument('ion_input'); parser.add_argument('ionization_input')
    args = parser.parse_args()
    report = run(args.destination, args.baseline, args.neutral_input, args.ion_input, args.ionization_input)
    raise SystemExit(0 if report['success'] else 1)
