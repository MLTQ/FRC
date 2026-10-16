"""Complete joint atomic-source evolutions and existing independent physical limits."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV
from plasma_bench.adf07_ionization import read_ionization
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.collisional_levels import CollisionalLevels
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from .impact_ionization import hydrogen_reference
from .radiative_recombination import power_law_reference
from .implicit_plasma_records import provenance, freeze


def evolve(model, state, density, duration, steps, deposit=0., pulse=1e-6):
    times = np.linspace(0, duration, steps+1); states = [state.copy()]; failure = None
    def heat(t):
        return 2*deposit/pulse*np.sin(np.pi*t/pulse)**2 if 0 < t < pulse else 0.
    for step, (start, end) in enumerate(zip(times[:-1], times[1:]), 1):
        try:
            candidate = model.step(state, end-start, density,
                initial_heat=heat(start), midpoint_heat=heat(.5*(start+end)))
        except (ValueError, RuntimeError) as error:
            failure = dict(step=step, kind=type(error).__name__, message=str(error)); break
        state = candidate; states.append(state.copy())
    states = np.asarray(states); times = times[:len(states)]; populations = states[:, :model.count]
    number = states[:, model.electrons]; charge = populations@model.charge
    temperature = states[:, model.thermal]/(1.5*number)
    atomic = populations@model.energy_ev
    total = atomic+states[:, model.thermal]+states[:, model.photons]
    inside = np.minimum(times/pulse, 1.)
    analytic_heat = deposit*(inside-np.sin(2*np.pi*inside)/(2*np.pi))
    diagnostics = [model.generator(s, density)[1] for s in states]
    excluded = np.array([d['missing_excited_ionization_population'] for d in diagnostics])
    ionization = np.array([np.sum(d['ionization_events_per_nucleus_s']) for d in diagnostics])
    excited_ionization = np.array([d['resolved_excited_ionization_events_per_nucleus_s'] for d in diagnostics])
    rr_events = np.array([np.sum(d['recombination_events_per_nucleus_s']) for d in diagnostics])
    total_events = np.trapezoid(ionization, times)
    summary = dict(success=failure is None, failure=failure, completed_steps=len(states)-1,
        requested_steps=steps, dt_s=duration/steps, final_time_s=float(times[-1]),
        final_mean_charge=float(charge[-1]), final_electron_number=float(number[-1]),
        final_temperature_ev=float(temperature[-1]), final_photon_energy_ev=float(states[-1, model.photons]),
        final_excitation_energy_ev=float(atomic[-1]-populations[-1]@(model.ionization.atomic_energy[model.charge]/EV)),
        maximum_number_error=float(abs(populations.sum(axis=1)-populations[0].sum()).max()),
        maximum_charge_balance_error=float(abs(number-charge).max()),
        maximum_energy_balance_ev=float(abs(total-total[0]-states[:, model.heating]).max()),
        maximum_energy_minus_analytic_heat_ev=float(abs(total-total[0]-analytic_heat).max()),
        minimum_population=float(populations.min()),
        maximum_missing_excited_ionization_population=float(excluded.max()),
        maximum_missing_excited_parent_rr_population=max(d['missing_excited_parent_rr_population'] for d in diagnostics),
        integrated_ionization_events_per_nucleus=float(total_events),
        integrated_excited_ionization_event_fraction=float(np.trapezoid(excited_ionization,times)/total_events) if total_events else 0.,
        integrated_recombination_events_per_nucleus=float(np.trapezoid(rr_events,times)),
        maximum_outside_ionization_fit_event_fraction=max(d['outside_ionization_fit_event_fraction'] for d in diagnostics),
        maximum_discrepant_threshold_event_fraction=max(d['ionization_threshold_discrepancy_over_5pct_event_fraction'] for d in diagnostics),
        maximum_unknown_rr_range_event_fraction=max(d['unknown_rr_range_event_fraction'] for d in diagnostics),
        maximum_outside_known_rr_range_event_fraction=max(d['outside_known_rr_range_event_fraction'] for d in diagnostics),
        maximum_ionization_threshold_discrepancy_ev=max(d['ionization_threshold_discrepancy_max_ev'] for d in diagnostics))
    return summary, dict(time_s=times, states=states, mean_charge=charge, temperature_ev=temperature,
        total_energy_ev=total, analytic_heat_ev=analytic_heat,
        energy_balance_ev=total-total[0]-states[:, model.heating],
        missing_excited_ionization_population=excluded, ionization_events_per_nucleus_s=ionization,
        excited_ionization_events_per_nucleus_s=excited_ionization)


def argon_models(neutral_path, ion_path, ionization_path):
    data = Path(__file__).resolve().parents[1]/'data/atomic'
    ionization = ImpactIonization.from_file(data/'argon/argon.json')
    recombination = RadiativeRecombination.from_file(data/'recombination/argon_rr.json')
    neutral = CollisionalLevels(neutral_path,
        radiative_path=data/'excitation/argon_low_radiative.json', exclude_processes=('S',))
    ion = CollisionalLevels(ion_path); resolved = read_ionization(ionization_path)
    mappings = []; identities = []
    for metastable, configuration, multiplicity, orbital, weight in (
            (1,'3p5(2p)',2,1,6), (2,'3p4(3p)3d(4d)',4,2,20), (3,'3p4(3p)3d(4f)',4,3,28)):
        matches = [k for k, level in enumerate(ion.data['levels'])
                   if level['configuration'].lower()==configuration and level['multiplicity']==multiplicity
                   and level['orbital_L']==orbital and level['degeneracy']==weight]
        rows = [row for row in resolved['rows'] if row['initial_metastable']==metastable]
        if len(matches)!=1 or len(rows)!=1:
            raise ValueError('Unique inspected Ar II LS-term mapping required')
        mappings.append((1,matches[0],rows[0]))
        identities.append(dict(source_metastable=metastable, level_index=matches[0],
            configuration=configuration, source_effective_potential_cm1=rows[0]['effective_potential_cm1']))
    models = {'joint': AtomicNetwork(ionization,recombination,{0:neutral,1:ion},resolved_ionization=mappings),
              'ground_only': AtomicNetwork(ionization,recombination),
              'joint_without_metastable_ionization': AtomicNetwork(ionization,recombination,
                    {0:neutral,1:ion},resolved_ionization=mappings[:1])}
    inputs = dict(neutral=dict(path=neutral.data['input_path'],sha256=neutral.data['input_sha256']),
        ion=dict(path=ion.data['input_path'],sha256=ion.data['input_sha256']),
        metastable_ionization=dict(path=resolved['input_path'],sha256=resolved['input_sha256'], mappings=identities),
        ground_ionization_sha256=ionization.input_sha256, ground_rr_sha256=recombination.input_sha256,
        neutral_radiation_sha256=neutral.radiation.input_sha256)
    return models, inputs


def run(destination, neutral_path, ion_path, ionization_path):
    out = Path(destination)
    if out.exists(): raise ValueError('Fresh output directory required')
    models, inputs = argon_models(neutral_path, ion_path, ionization_path)
    out.mkdir(parents=True); started = time.monotonic()
    identity = provenance(['examples.atomic_network']); freeze(out, identity)
    report = dict(terminal=False,source_sha256=identity['source_sha256'],inputs=inputs,cases=[],references=[],comparisons=[],
        scope='Homogeneous joint level/charge/electron-number/energy evolution. Supplied ground and selected Ar II metastable ionization, ground-parent RR with immediate ground cascade, resolved excitation and separate photon production. Missing excited channels, unresolved ionization products, absorption, transfer, neutrals, momentum and moving coupling preclude a full plasma cooling or formation claim.',
        input_policy='Personal OPEN-ADAS inputs outside repository; no raw tables or reconstructed rate arrays archived.')
    histories = {}
    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    # Reuse independently derived existing analytic references, not solver internals.
    hydrogen = ImpactIonization(dict(nuclear_charge=1,transitions=[dict(charge=0,
        threshold_fit_ev=13.6,ionization_ev=13.6,A_cm3_s=2.91e-8,P=0.,X=.232,K=.39,
        minimum_temperature_ev=1.,maximum_temperature_ev=20000.)]))
    power_rr = RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,
        kind='power_law',A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=.7,
        minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    for name, duration, steps_list, fraction, temperature in (
            ('analytic_di',20e-6,(100,200,400),[.8,.2],30.),
            ('analytic_rr',.03,(50,100,200),[.2,.8],2.)):
        model = AtomicNetwork(hydrogen,power_rr,ionization_enabled=name=='analytic_di',
                              recombination_enabled=name=='analytic_rr')
        for steps in steps_list:
            summary, history = evolve(model,model.initial(fraction,temperature),1e20,duration,steps)
            if name=='analytic_di':
                charge, te = hydrogen_reference(history['time_s'],1e20,.2,30.)
            else:
                charge, te, photon = power_law_reference(history['time_s'])
                history['reference_photon_ev'] = photon
                summary['maximum_photon_error_ev'] = float(abs(history['states'][:,model.photons]-photon).max())
            history.update(reference_charge=charge,reference_temperature_ev=te)
            summary.update(name=name,steps=steps,maximum_reference_charge_error=float(abs(history['mean_charge']-charge).max()),
                maximum_reference_temperature_error_ev=float(abs(history['temperature_ev']-te).max()))
            report['references'].append(summary); np.savez_compressed(out/f'{name}_{steps}.npz',**history);save()
            print(name,steps,summary['success'],summary['maximum_reference_charge_error'],flush=True)
    # Fully neutral radiative decay is a physical zero-electron limit, not a seed-floor case.
    model=models['joint']; radiation=model.levels[0].radiation; local_upper=2
    outgoing=np.flatnonzero(radiation.upper==local_upper)
    if len(outgoing)!=1 or radiation.lower[outgoing[0]]!=0:
        raise ValueError('Independent neutral exponential reference requires one decay to ground')
    rate=radiation.A[outgoing[0]]; upper=model.blocks[0][local_upper]
    fractions=np.zeros(19); fractions[0]=1.; state=model.initial(fractions,0.)
    state[model.ground[0]]=0.; state[upper]=1.; states=[state.copy()]
    times=np.linspace(0,40e-9,101)
    for start,end in zip(times[:-1],times[1:]):
        state=model.step(state,end-start,1e20);states.append(state.copy())
    states=np.asarray(states); exact=np.exp(-rate*times)
    spacing=model.energy_ev[upper]-model.energy_ev[model.ground[0]]
    report['references'].append(dict(name='neutral_decay_zero_electrons',steps=100,success=True,
        maximum_population_error=float(abs(states[:,upper]-exact).max()),
        maximum_photon_error_ev=float(abs(states[:,model.photons]-spacing*(1-exact)).max()),
        maximum_electron_number=float(abs(states[:,model.electrons]).max()),
        maximum_electron_energy_ev=float(abs(states[:,model.thermal]).max())))
    np.savez_compressed(out/'neutral_decay_zero_electrons.npz',time_s=times,states=states,
                        reference_upper_population=exact,reference_photon_ev=spacing*(1-exact));save()
    density = 8.961429541213399e19
    for case,seed,deposit in [('singly_ionized',1.,0.),('partial_seed',.1,0.),('heated_seed',.1,20.)]:
        fraction=np.zeros(19);fraction[0]=1-seed;fraction[1]=seed
        for name,model in models.items():
            pair=[]
            for steps in ((700,1400) if name=='joint' else (1400,)):
                summary,history=evolve(model,model.initial(fraction,10.),density,3.5e-6,steps,deposit)
                key=f'{case}_{name}_{steps}';histories[key]=history
                np.savez_compressed(out/f'{key}.npz',**history)
                summary.update(case=case,model=name,steps=steps,nucleus_density_m3=density,
                    initial_singly_ionized_fraction=seed,initial_temperature_ev=10.,deposit_ev_per_nucleus=deposit)
                report['cases'].append(summary);pair.append((summary,history));save()
                print(case,name,steps,summary['success'],summary['final_mean_charge'],summary['final_temperature_ev'],flush=True)
            if len(pair)==2 and all(r[0]['success'] for r in pair):
                a,b=pair[0][1],pair[1][1]
                report['comparisons'].append(dict(case=case,comparison='joint_time_halving_5_to_2.5ns',
                    maximum_population_l1_difference=float(np.sum(abs(a['states'][:,:model.count]-b['states'][::2,:model.count]),axis=1).max()),
                    maximum_temperature_difference_ev=float(abs(a['temperature_ev']-b['temperature_ev'][::2]).max()),
                    maximum_charge_difference=float(abs(a['mean_charge']-b['mean_charge'][::2]).max())))
                save()
    fig,axes=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
    for column,case in enumerate(('singly_ionized','partial_seed','heated_seed')):
        for name in models:
            data=histories[f'{case}_{name}_1400'];t=data['time_s']*1e6
            label=name.replace('_',' ')
            axes[0,column].plot(t,data['temperature_ev'],label=label)
            axes[1,column].plot(t,data['mean_charge'],label=label)
            axes[2,column].plot(t,data['missing_excited_ionization_population'],label=label)
        axes[0,column].set_title(case.replace('_',' '))
        for row,label in enumerate(('Electron temperature [eV]','Electrons per nucleus','Population lacking excited ionization rates')):
            axes[row,column].set(xlabel='Time [µs]',ylabel=label);axes[row,column].grid(alpha=.2)
            axes[row,column].legend(fontsize=7,frameon=False)
    fig.suptitle('Joint charge, level and electron-energy evolution\nIncomplete atomic channels and unresolved products remain explicit; photons are produced, not established escape')
    fig.savefig(out/'joint_atomic_evolution.png',dpi=170);fig.savefig(out/'joint_atomic_evolution.pdf');plt.close(fig)
    report['terminal']=True;save();print('Joint atomic source complete',out,flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('neutral_input')
    parser.add_argument('ion_input');parser.add_argument('ionization_input');args=parser.parse_args()
    run(args.destination,args.neutral_input,args.ion_input,args.ionization_input)
