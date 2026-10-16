"""Full physical source histories for array-based moment assembly."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization, evolve, conversion_reference
from .reactive_material import recombination
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh source comparison output required')
    identity = provenance(['examples.atomic_vector_moments'])
    out.mkdir(parents=True)
    freeze(out, identity)
    began = time.perf_counter()
    report = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[], comparisons=[],
        scope='Same finite-mass reaction physics, exponential midpoint and complete durations. '
              'Independent constant-rate conversion limits plus actual ground-ladder argon source; '
              'no density coarsening, moving-field or formation claim.')
    def save():
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    def trajectory(name, model, initial, density, duration, steps, deposit=0., analytic=None):
        started = time.perf_counter()
        row, history = evolve(model, initial, density, duration, steps, deposit)
        row.update(name=name, elapsed_seconds=time.perf_counter()-started,
                   completed_time_s=float(history['time_s'][-1]))
        if analytic is not None:
            exact = analytic(history['time_s'])
            history['analytic_states'] = exact
            delta = history['states']-exact
            row.update(maximum_analytic_population_error=float(abs(delta[:, :model.network.count]).max()),
                       maximum_analytic_scaled_momentum_error=float(abs(delta[:, model.momentum]).max()),
                       maximum_analytic_translational_energy_error_ev=float(abs(delta[:, model.energy]).max()))
        np.savez_compressed(out/(name+'.npz'), **history)
        report['cases'].append(row)
        save()
        print(name, 'success', row['success'], 'seconds', row['elapsed_seconds'],
              'failure', row['failure'], flush=True)
        return row, history
    save()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for capture in (False, True):
        network = AtomicNetwork(ConstantIonization(), recombination(),
                                ionization_enabled=not capture, recombination_enabled=capture, radiation=False)
        model = VectorizedAtomicMoments(network, 1.00794*atomic_mass, [0., 0., 0.])
        fraction, temperature = (.8, 2.) if capture else (.2, 30.)
        initial = model.initial(network.initial([1-fraction, fraction], temperature), .1,
            np.array([[1000., 2000., -500.], [-1000., 500., 0.], [5000., -3000., 500.]]))
        errors = []
        for steps in (50, 100, 200):
            row, history = trajectory(f"{'rr' if capture else 'di'}_{steps}", model, initial,
                1e20, .02 if capture else 2e-6, steps,
                analytic=lambda t: conversion_reference(model, initial, t, 1e20, capture))
            errors.append(row['maximum_analytic_population_error'])
        axes[0].loglog([50, 100, 200], errors, 'o-', label='Recombination' if capture else 'Ionization')

    ion = ImpactIonization.from_file('data/atomic/argon/argon.json')
    rr = RadiativeRecombination.from_file('data/atomic/recombination/argon_rr.json')
    network = AtomicNetwork(ion, rr)
    report['inputs'] = dict(ionization_sha256=ion.input_sha256, recombination_sha256=rr.input_sha256)
    velocity = np.zeros((network.count+1, 3))
    velocity[:-1, 2] = np.where(network.charge == 0, 1000., -1000.)
    velocity[-1] = [0., 10000., 0.]
    atomic = network.initial(np.r_[.9, .1, np.zeros(17)], 10.)
    density = 8.961429541213399e19
    histories = {}
    for steps in (700, 1400):
        for name, constructor in [('scalar', FixedFrameAtomicMoments), ('vector', VectorizedAtomicMoments)]:
            model = constructor(network, 39.948*atomic_mass, [0., 0., 0.])
            initial = model.initial(atomic, .1, velocity)
            row, history = trajectory(f'argon_{name}_{steps}', model, initial, density, 3.5e-6, steps, 20.)
            histories[name, steps] = (row, history)
            axes[1].plot(history['time_s']*1e6, history['electron_temperature_ev'],
                         '--' if name == 'scalar' else '-', label=f'{name}, {steps}')
            axes[2].plot(history['time_s']*1e6,
                         history['energy_budget_ev']-history['energy_budget_ev'][0], label=f'{name}, {steps}')
        old, new = histories['scalar', steps], histories['vector', steps]
        if old[0]['success'] and new[0]['success']:
            delta = new[1]['states']-old[1]['states']
            report['comparisons'].append(dict(steps=steps,
                trajectory_speedup=old[0]['elapsed_seconds']/new[0]['elapsed_seconds'],
                maximum_population_difference=float(abs(delta[:, :network.count]).max()),
                maximum_momentum_difference_scaled=float(abs(delta[:, model.momentum]).max()),
                maximum_energy_difference_ev=float(abs(delta[:, model.energy]).max()),
                maximum_current_difference_a_m2=float(abs(new[1]['current_a_m2']-old[1]['current_a_m2']).max())))
        save()
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    axes[0].set(xlabel='Complete-history time steps', ylabel='Population error vs analytic solution')
    axes[1].set(xlabel='Time [µs]', ylabel='Electron temperature [eV]')
    axes[2].set(xlabel='Time [µs]', ylabel='Total-energy error [eV / initial nucleus]')
    fig.savefig(out/'reaction_assembly.png', dpi=160)
    fig.savefig(out/'reaction_assembly.pdf')
    plt.close(fig)
    if provenance(['examples.atomic_vector_moments']) != identity:
        raise RuntimeError('Source changed during complete reaction histories')
    report['terminal'] = True
    report['success'] = all(c['success'] for c in report['cases'])
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    args = parser.parse_args()
    result = run(args.destination)
    raise SystemExit(0 if result['success'] else 1)
