"""Complete reactive momentum/energy trajectories and independent conversion limits."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV, electron_mass as ME
from plasma_bench.atomic_moments import AtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_network import argon_models
from .implicit_plasma_records import provenance, freeze


class ConstantIonization:
    """Declared mathematical constant-coefficient limit, not an atomic fit."""
    nuclear_charge = 1
    charge = np.arange(2)
    atomic_energy = np.array([0., 13.6])*EV
    fit_energy = np.array([13.6])*EV

    def coefficients(self, temperature):
        return np.array([1e-14]), np.array([True])


def evolve(model, initial, density, duration, steps, deposit=0.):
    state = initial.copy(); times = np.linspace(0, duration, steps+1)
    states = [state.copy()]; failure = None; drift = []; thermal_ratio = []
    def heat(t):return 2*deposit/1e-6*np.sin(np.pi*t/1e-6)**2 if 0 < t < 1e-6 else 0.
    for index, (start, end) in enumerate(zip(times[:-1], times[1:]), 1):
        try:
            candidate = model.step(state, end-start, density,
                initial_heat=heat(start), midpoint_heat=heat(.5*(start+end)))
        except (ValueError, RuntimeError) as error:
            failure = dict(step=index, kind=type(error).__name__, message=str(error)); break
        state = candidate; states.append(state.copy())
    states = np.asarray(states); times = times[:len(states)]
    rows = [model.observe(s, density) for s in states]
    history = {key:np.asarray([row[key] for row in rows]) for key in rows[0]}
    for s in states:
        _, d = model.generator(s, density)
        drift.append(d['maximum_drift_over_electron_thermal'])
        thermal_ratio.append(d['maximum_heavy_thermal_over_electron_thermal'])
    summary = dict(success=failure is None, failure=failure, steps=steps,
        accepted_steps=len(states)-1, duration_s=duration, dt_s=duration/steps,
        maximum_number_error=float(abs(history['number']-history['number'][0]).max()),
        maximum_charge_error=float(abs(history['charge_balance']).max()),
        maximum_mass_error=float(abs(history['mass_over_neutral']-history['mass_over_neutral'][0]).max()),
        maximum_scaled_momentum_error=float(np.linalg.norm(history['momentum_scaled']-history['momentum_scaled'][0],axis=1).max()),
        maximum_energy_error_ev=float(abs(history['energy_budget_ev']-history['energy_budget_ev'][0]).max()),
        maximum_drift_over_electron_thermal=max((x for x in drift if x is not None),default=None),
        maximum_heavy_thermal_over_electron_thermal=max((x for x in thermal_ratio if x is not None),default=None),
        final_electron_number=float(history['electrons'][-1]),
        final_electron_temperature_ev=float(history['electron_temperature_ev'][-1]),
        final_current_a_m2=history['current_a_m2'][-1].tolist(),
        final_photon_energy_ev=float(history['photon_energy_ev'][-1]))
    return summary, dict(time_s=times, states=states, **history)


def conversion_reference(model, initial, times, density, recombination):
    """Closed logistic/reciprocal number law and per-event incoming/outgoing moments."""
    n = model.network; f = model.fields(initial); x0 = initial[n.electrons]
    if recombination:
        x = x0/(1+density*1e-18*x0*times); converted = x0-x
    else:
        x = 1/(1+(1/x0-1)*np.exp(-density*1e-14*times)); converted = x-x0
    expected = np.repeat(initial[None, :], len(times), axis=0)
    expected[:, :2] = np.column_stack((1-x,x)); expected[:, n.electrons] = x
    masses = model.mass; velocities = f['velocity']; energies = initial[model.energy]
    if not recombination:
        parent_energy = energies[0]/(1-x0)
        expected[:, model.momentum[0]] = initial[model.momentum[0]]*(1-x)[:,None]/(1-x0)
        expected[:, model.heavy_energy[0]] = energies[0]*(1-x)/(1-x0)
        for child, fraction in ((1,masses[1]/masses[0]),(2,ME/masses[0])):
            expected[:, model.momentum[child]] += converted[:,None]*masses[child]*velocities[0]/model.momentum_scale
            expected[:, model.energy[child]] += converted*fraction*parent_energy
        expected[:, n.thermal] -= 13.6*converted
    else:
        vi,ve = velocities[1:]; mi,mn = masses[1],masses[0]
        ui,ue = f['thermal_ev'][1:]/x0
        pair_momentum = mi*vi+ME*ve
        product = .5*(pair_momentum@pair_momentum)/(mn*EV)+mi/mn*ui+ME/mn*ue
        photons = 13.6+.5*mi*ME/mn*np.sum((vi-ve)**2)/EV+ME/mn*ui+mi/mn*ue
        for child in (1,2):
            expected[:, model.momentum[child]] = initial[model.momentum[child]]*x[:,None]/x0
            expected[:, model.energy[child]] = energies[child]*x/x0
        expected[:, model.momentum[0]] += converted[:,None]*pair_momentum/model.momentum_scale
        expected[:, model.heavy_energy[0]] += converted*product
        expected[:, n.photons] += converted*photons
    return expected


def run(destination, neutral_path, ion_path, ionization_path):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh output directory required')
    models,inputs=argon_models(neutral_path,ion_path,ionization_path)
    identity=provenance(['examples.atomic_moments']);out.mkdir(parents=True);freeze(out,identity)
    started=time.monotonic();histories={};report=dict(terminal=False,source_sha256=identity['source_sha256'],
        inputs=inputs,cases=[],comparisons=[],scope='Homogeneous reacting level/electron fluids with finite masses and momentum/energy transfer. Small relative drift, isotropic capture and leading nonrelativistic photon recoil approximation. No electromagnetic/spatial coupling or photon escape claim. Original missing atomic channels and range limits remain.',
        input_policy='Private provider tables and reconstructed rate matrices are not copied into outputs.')
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def record(name,model,initial,density,duration,steps,deposit=0.,reference=None):
        summary,history=evolve(model,initial,density,duration,steps,deposit)
        summary['name']=name
        if reference is not None:
            exact=reference(history['time_s']); history['analytic_states']=exact
            delta=history['states']-exact;n=model.network
            summary['maximum_analytic_population_error']=float(abs(delta[:,:n.count]).max())
            summary['maximum_analytic_scaled_momentum_error']=float(abs(delta[:,model.momentum]).max())
            summary['maximum_analytic_translational_energy_error_ev']=float(abs(delta[:,model.energy]).max())
        np.savez_compressed(out/f'{name}.npz',**history)
        report['cases'].append(summary);histories[name]=(model,history);save()
        print('Atomic moments',name,summary['success'],summary['maximum_energy_error_ev'],flush=True)
    save()
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    for recombination in (False,True):
        network=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=not recombination,
                              recombination_enabled=recombination,radiation=False)
        model=AtomicMoments(network,1.00794*atomic_mass)
        x0=.8 if recombination else .2;te=2. if recombination else 30.
        initial=model.initial(network.initial([1-x0,x0],te),.1,
                              np.array([[1000.,2000.,-500.],[-1000.,500.,0.],[5000.,-3000.,500.]]))
        duration=.02 if recombination else 2e-6
        for steps in (50,100,200):
            record(f"analytic_{'rr' if recombination else 'di'}_{steps}",model,initial,1e20,duration,steps,
                   reference=lambda times:conversion_reference(model,initial,times,1e20,recombination))
    network=models['joint'];model=AtomicMoments(network,39.948*atomic_mass);density=8.961429541213399e19
    velocity=np.zeros((network.count+1,3));velocity[:-1,2]=np.where(network.charge==0,1000.,-1000.)
    velocity[-1]=[0.,10000.,0.]
    for name,seed,deposit in [('all_ion',1.,0.),('seed',.1,0.),('heated_seed',.1,20.)]:
        fractions=np.zeros(19);fractions[:2]=[1-seed,seed]
        atomic=network.initial(fractions,10.)
        initial=model.initial(atomic,.1,velocity)
        for steps in (700,1400):record(f'{name}_{steps}',model,initial,density,3.5e-6,steps,deposit)
        if name=='heated_seed':
            boosted=model.initial(atomic,.1,velocity+np.array([20000.,-10000.,5000.]))
            record('heated_seed_boosted_1400',model,boosted,density,3.5e-6,1400,deposit)
    # Moving excited neutral decay: exactly no free electrons or electron energy.
    neutral_network=AtomicNetwork(network.ionization,network.recombination,{0:network.levels[0]},
                                  ionization_enabled=False,recombination_enabled=False)
    neutral_model=AtomicMoments(neutral_network,39.948*atomic_mass)
    level=neutral_network.levels[0];upper=int(level.radiation.upper[0]);lower=int(level.radiation.lower[0])
    outgoing=np.flatnonzero(level.radiation.upper==upper)
    if len(outgoing)!=1:raise ValueError('Independent single-line decay requires one radiative destination')
    atomic=neutral_network.initial(np.r_[1.,np.zeros(18)],0.)
    atomic[neutral_network.ground[0]]=0.;atomic[upper]=1.
    neutral_initial=neutral_model.initial(atomic,.02585,np.tile([1000.,-2000.,500.],(neutral_network.count+1,1)))
    def decay_reference(times):
        fraction=np.exp(-level.radiation.A[outgoing[0]]*times)
        exact=np.repeat(neutral_initial[None,:],len(times),axis=0)
        exact[:,upper]=fraction;exact[:,lower]=1-fraction
        for indices in (neutral_model.momentum,neutral_model.heavy_energy):
            original=neutral_initial[indices[upper]]
            shape=(-1,1) if np.ndim(original) else (-1,)
            exact[:,indices[upper]]=fraction.reshape(shape)*original
            exact[:,indices[lower]]=(1-fraction).reshape(shape)*original
        exact[:,neutral_network.photons]=(1-fraction)*(neutral_network.energy_ev[upper]-neutral_network.energy_ev[lower])
        return exact
    record('moving_neutral_decay',neutral_model,neutral_initial,density,40e-9,80,reference=decay_reference)
    for name in ('all_ion','seed','heated_seed'):
        _,coarse=histories[f'{name}_700'];_,fine=histories[f'{name}_1400']
        if len(coarse['states'])!=701 or len(fine['states'])!=1401:
            report['comparisons'].append(dict(name=name,available=False,reason='A complete trajectory failed; no comparison of unequal prefixes.'))
            continue
        report['comparisons'].append(dict(name=name,
            maximum_population_l1=float(abs(coarse['states'][:,:network.count]-fine['states'][::2,:network.count]).sum(axis=1).max()),
            maximum_temperature_difference_ev=float(abs(coarse['electron_temperature_ev']-fine['electron_temperature_ev'][::2]).max()),
            maximum_current_difference_a_m2=float(np.linalg.norm(coarse['current_a_m2']-fine['current_a_m2'][::2],axis=1).max())))
    _,base=histories['heated_seed_1400'];_,boost=histories['heated_seed_boosted_1400']
    if len(base['states'])==len(boost['states'])==1401:
        report['galilean_comparison']=dict(
            maximum_population_l1=float(abs(base['states'][:,:network.count]-boost['states'][:,:network.count]).sum(axis=1).max()),
            maximum_temperature_difference_ev=float(abs(base['electron_temperature_ev']-boost['electron_temperature_ev']).max()),
            maximum_current_difference_a_m2=float(np.linalg.norm(base['current_a_m2']-boost['current_a_m2'],axis=1).max()))
    else:
        report['galilean_comparison']=dict(available=False,reason='A complete trajectory failed; no comparison of unequal prefixes.')
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for name in ('all_ion','seed','heated_seed'):
        _,a=histories[f'{name}_1400'];t=a['time_s']*1e6
        for ax,key in zip(axes.flat,('electrons','electron_temperature_ev','current_a_m2','energy_budget_ev')):
            value=a[key][:,1] if key=='current_a_m2' else a[key]-a[key][0] if key=='energy_budget_ev' else a[key]
            ax.plot(t,value,label=name+(' (incomplete)' if len(t)!=1401 else ''))
    for ax,label in zip(axes.flat,('Electrons / nucleus','Electron temperature [eV]','Azimuthal current [A/m²]','Total energy budget residual [eV/nucleus]')):
        ax.set(xlabel='Time [µs]',ylabel=label);ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Homogeneous reactive moments: small-drift closure; no field evolution')
    fig.savefig(out/'reactive_moments.png',dpi=170);fig.savefig(out/'reactive_moments.pdf');plt.close(fig)
    if provenance(['examples.atomic_moments'])!=identity:raise RuntimeError('Source changed during evolution')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save()
    print('Atomic moments terminal',report['success'],out,flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('neutral_input')
    parser.add_argument('ion_input');parser.add_argument('ionization_input');args=parser.parse_args()
    report=run(args.destination,args.neutral_input,args.ion_input,args.ionization_input)
    raise SystemExit(0 if report['success'] else 1)
