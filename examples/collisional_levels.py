"""Complete excitation/thermal evolutions with supplied personal argon inputs."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E, Boltzmann as KB
from scipy.optimize import brentq
from plasma_bench.collisional_levels import CollisionalLevels
from .implicit_plasma_records import provenance, freeze


def boltzmann(model, temperature):
    weight = model.weight*np.exp(-model.energy_ev/temperature)
    return weight/weight.sum()


def evolve(model, initial_temperature, duration, steps, radiation):
    count = len(model.energy_ev); electrons = 1.; density = 1e20
    state = np.zeros(count+2); state[0] = 1.; state[count] = 1.5*electrons*initial_temperature
    initial_energy = state[count]; states = [state.copy()]; failure = None
    for step in range(1, steps+1):
        try:
            state = model.step(state, duration/steps, density, electrons, radiation=radiation)
        except (ValueError, RuntimeError) as error:
            failure = dict(step=step, type=type(error).__name__, message=str(error))
            break
        states.append(state.copy())
    states = np.array(states); populations = states[:, :count]
    atomic = populations@model.energy_ev; temperature = states[:, count]/(1.5*electrons)
    total = atomic+states[:, count]+states[:, count+1]
    atomic_entropy = -np.sum(np.where(populations > 0, populations*np.log(np.where(populations > 0, populations/model.weight, 1.)), 0.), axis=1)
    entropy = atomic_entropy+1.5*electrons*np.log(temperature)
    summary = dict(success=failure is None, failure=failure, requested_steps=steps,
        completed_steps=len(states)-1, duration_s=duration, dt_s=duration/steps,
        initial_temperature_ev=initial_temperature, electron_density_m3=density,
        electrons_per_modeled_atom=electrons, spontaneous_radiation=radiation,
        final_temperature_ev=float(temperature[-1]), final_excitation_ev_per_atom=float(atomic[-1]),
        produced_photon_ev_per_atom=float(states[-1, count+1]),
        max_population_sum_error=float(np.max(abs(populations.sum(axis=1)-1))),
        max_relative_energy_change=float(np.max(abs(total/initial_energy-1))),
        minimum_population=float(populations.min()),
        minimum_sampled_matter_entropy_increment_over_kB=float(np.diff(entropy).min()) if len(entropy)>1 else None)
    if not radiation:
        low = model.temperature[0]*KB/E
        equilibrium_temperature = brentq(lambda t: 1.5*electrons*t+boltzmann(model,t)@model.energy_ev-initial_energy,
                                         low, initial_temperature, xtol=1e-13)
        equilibrium = boltzmann(model, equilibrium_temperature)
        summary.update(boltzmann_temperature_ev=equilibrium_temperature,
            final_temperature_error_ev=float(temperature[-1]-equilibrium_temperature),
            final_population_l1_error=float(np.sum(abs(populations[-1]-equilibrium))),
            boltzmann_populations=equilibrium.tolist())
    return summary, dict(time_s=np.arange(len(states))*duration/steps, states=states,
                        temperature_ev=temperature, excitation_ev=atomic, total_energy_ev=total,
                        matter_entropy_over_kB=entropy)


def run(destination, neutral_path, ion_path):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh output directory required')
    identity = provenance(['examples.collisional_levels']); started = time.monotonic()
    models = {'ArI': CollisionalLevels(neutral_path,
        radiative_path='data/atomic/excitation/argon_low_radiative.json', exclude_processes=('S',)),
        'ArII_terms': CollisionalLevels(ion_path)}
    out.mkdir(parents=True); freeze(out, identity)
    record = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[], comparisons=[],
        scope='Homogeneous fixed-charge Maxwellian excitation component with finite electron heat and separate produced-photon energy. No ionization, absorption, moving plasma, escape or validated net cooling closure. ArI electrons are balanced by an unevolved background ion population; ArII is singly ionized with ne/nion=1.',
        input_policy='Explicit personal local OPEN-ADAS files outside the repository; no raw tables or reconstructed rate arrays are copied into these outputs.',
        inputs={name:dict(path=m.data['input_path'],sha256=m.data['input_sha256'],levels=len(m.energy_ev),
            excitation_pairs=len(m.lower), temperature_range_k=m.temperature[[0,-1]].tolist(),
            radiative_input_sha256=m.radiation.input_sha256, excluded_processes=m.excluded_processes)
            for name,m in models.items()})
    histories = {}
    def save():
        record['elapsed_seconds'] = time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    save()
    for name, model in models.items():
        for radiation in (True, False):
            temperatures = (5.,10.,20.) if radiation else (10.,)
            duration = 3.5e-6 if radiation else 1e-3
            for temperature in temperatures:
                pair = []
                for steps in ((350,700) if radiation else (10000,20000)):
                    key = f'{name}_{"radiation" if radiation else "closed"}_{temperature:g}ev_{steps}'
                    summary, data = evolve(model,temperature,duration,steps,radiation)
                    summary.update(name=key, ion=name); record['cases'].append(summary)
                    np.savez_compressed(out/f'{key}.npz', **data)
                    histories[key] = data; pair.append((summary,data)); save()
                    print('Resolved excitation',key,summary['success'],summary['completed_steps'],flush=True)
                if all(item[0]['success'] for item in pair):
                    first, second = [item[1] for item in pair]
                    record['comparisons'].append(dict(ion=name, radiation=radiation, initial_temperature_ev=temperature,
                        max_population_l1_difference=float(np.max(np.sum(abs(first['states'][:,:len(model.energy_ev)]-second['states'][::2,:len(model.energy_ev)]),axis=1))),
                        max_temperature_difference_ev=float(np.max(abs(first['temperature_ev']-second['temperature_ev'][::2])))))
                    save()
    fig, axes = plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    equilibrium_fig, eq_axes = plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for column, (name,model) in enumerate(models.items()):
        for temperature in (5.,10.,20.):
            data=histories[f'{name}_radiation_{temperature:g}ev_700']
            axes[0,column].plot(data['time_s']*1e6,data['temperature_ev'],label=f'Initially {temperature:g} eV')
        data=histories[f'{name}_radiation_10ev_700'];count=len(model.energy_ev)
        for values,label in ((data['states'][:,count],'Electron thermal'),(data['excitation_ev'],'Atomic excitation'),
                             (data['states'][:,count+1],'Produced photons'),(data['total_energy_ev'],'Total')):
            axes[1,column].plot(data['time_s']*1e6,values,label=label)
        axes[0,column].set(title=name+' — radiation produced, no reabsorption',xlabel='Time [µs]',ylabel='Electron temperature [eV]')
        axes[1,column].set(xlabel='Time [µs]',ylabel='Energy [eV per modeled atom]')
        data=histories[f'{name}_closed_10ev_20000']
        result=next(r for r in record['cases'] if r['name']==f'{name}_closed_10ev_20000')
        eq_axes[0,column].plot(data['time_s']*1e6,data['temperature_ev'],label='Evolved electron temperature')
        eq_axes[0,column].axhline(result['boltzmann_temperature_ev'],color='k',ls='--',label='Independent energy + Boltzmann solution')
        eq_axes[0,column].set(title=name+' — closed collisional relaxation',xlabel='Time [µs]',ylabel='Electron temperature [eV]')
        eq_axes[1,column].plot(np.arange(count),data['states'][-1,:count],'o',label='Evolved populations')
        eq_axes[1,column].plot(np.arange(count),result['boltzmann_populations'],'x',label='Boltzmann populations')
        eq_axes[1,column].set(xlabel='Source level/term index',ylabel='Population per modeled atom')
    for ax in (*axes.flat,*eq_axes.flat):
        ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Resolved electron-impact excitation and finite electron energy\nFixed charge; produced radiation is not established escape or net plasma cooling')
    equilibrium_fig.suptitle('Independent finite-reservoir equilibrium limit\nSpontaneous radiation explicitly disabled in this closed reference problem')
    fig.savefig(out/'excitation_energy.png',dpi=170);equilibrium_fig.savefig(out/'closed_equilibrium.png',dpi=170)
    plt.close(fig);plt.close(equilibrium_fig)
    record['terminal']=True;save();print('Resolved excitation study complete',out,flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('neutral_input');parser.add_argument('ion_input')
    args=parser.parse_args();run(args.destination,args.neutral_input,args.ion_input)
