"""Complete recombining source trajectories and independent analytic evolution."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination, RecombiningPlasma
from .impact_ionization import evolve as ionization_only
from .implicit_plasma_records import provenance, freeze


def evolve(model, density, fraction, temperature_ev, duration, steps, deposit_ev=0., pulse=1e-6):
    f = np.asarray(fraction,float)
    y = np.r_[f,1.5*(model.charge@f)*temperature_ev*EV,0.]
    times = np.linspace(0.,duration,steps+1); states = [y.copy()]
    def rhs(t,state):
        heat = 2*deposit_ev*EV/pulse*np.sin(np.pi*t/pulse)**2 if 0<t<pulse else 0.
        df,de,dp,_ = model.source(state[:-2],state[-2],density,heating_j_per_nucleus_s=heat)
        return np.r_[df,de,dp]
    for t,end in zip(times[:-1],times[1:]):
        h=end-t
        k1=rhs(t,y); k2=rhs(t+h/2,y+h*k1/2)
        k3=rhs(t+h/2,y+h*k2/2); k4=rhs(end,y+h*k3)
        y=y+h*(k1+2*k2+2*k3+k4)/6
        model.source(y[:-2],y[-2],density)
        states.append(y.copy())
    states=np.asarray(states);populations=states[:,:-2];charge=populations@model.charge
    temperature=states[:,-2]/(1.5*charge*EV)
    atomic=populations@model.atomic_energy
    inside=np.minimum(times/pulse,1.)
    heating=deposit_ev*EV*(inside-np.sin(2*np.pi*inside)/(2*np.pi))
    total=atomic+states[:,-2]+states[:,-1]
    diagnostics=[model.source(s[:-2],s[-2],density)[3] for s in states]
    return dict(time_s=times,state=states,charge=charge,temperature_ev=temperature,
        atomic_energy_j_per_nucleus=atomic,analytic_heating_j_per_nucleus=heating,
        energy_balance_j_per_nucleus=total-total[0]-heating,
        unknown_rr_range_event_fraction=np.array([d['unknown_rr_range_event_fraction'] or 0. for d in diagnostics]),
        outside_known_rr_range_event_fraction=np.array([d['outside_known_rr_range_event_fraction'] or 0. for d in diagnostics]),
        outside_ionization_fit_event_fraction=np.array([d['outside_ionization_fit_event_fraction'] or 0. for d in diagnostics]))


def power_law_reference(times):
    # Analytic physical source limit; no production rates or time integrator called.
    # sigma(E) proportional to E^(-s-1/2) gives alpha proportional to T^-s.
    density=1e20;alpha0=1e-18;x0=.8;t0=2.;s=.7;chi=13.6
    beta=2*s*s/3
    x=x0*(1+(1+beta)*density*alpha0*x0*times)**(-1/(1+beta))
    temperature=t0*(x/x0)**(-2*s/3)
    photons=1.5*x0*t0+chi*x0-(1.5*x*temperature+chi*x)
    return x,temperature,photons


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('A new output directory is required')
    out.mkdir(parents=True);started=time.monotonic()
    identity=provenance(['examples.radiative_recombination']);freeze(out,identity)
    data=Path(__file__).resolve().parents[1]/'data/atomic'
    shutil.copytree(data/'argon',out/'ionization_inputs')
    shutil.copytree(data/'recombination',out/'recombination_inputs')
    argon=ImpactIonization.from_file(data/'argon/argon.json')
    rr=RadiativeRecombination.from_file(data/'recombination/argon_rr.json')
    model=RecombiningPlasma(argon,rr)
    report=dict(terminal=False,source_sha256=identity['source_sha256'],
        ionization_input_sha256=argon.input_sha256,recombination_input_sha256=rr.input_sha256,
        scope='Fixed-volume DI+RR only, Maxwellian ground-state parents, immediate ground-state cascade and photon escape. No DR, three-body, excitation, opacity, motion, transport or preparation prediction.',
        reference=[],argon=[],original_rate_comparison=[])
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    # Published author-calculated rates are independent input values, not our fitted function.
    originals=json.loads((data/'recombination/argon_rr_original.json').read_text())
    report['original_rate_missing_modern_initial_charges']=originals['missing_modern_initial_charges']
    for row in originals['records']:
        q=row['initial_charge'];samples=np.asarray(row['samples'])
        predicted=np.array([rr.coefficients(t*EV)[0][q-1] for t in samples[:,0]])
        relative=predicted/(samples[:,1]*1e-6)-1
        report['original_rate_comparison'].append(dict(initial_charge=q,points=len(samples),
            maximum_absolute_relative_difference=float(abs(relative).max())))
        axes[1,2].semilogx(samples[:,0],100*relative,label=f'Ar{q}+')
        np.savez_compressed(out/f'published_rates_q{q}.npz',original=samples,predicted_m3_s=predicted,relative_difference=relative)
    fixture_ion=ImpactIonization(dict(nuclear_charge=1,transitions=[dict(charge=0,
        threshold_fit_ev=13.6,ionization_ev=13.6,A_cm3_s=2.91e-8,P=0.,X=.232,K=.39,
        minimum_temperature_ev=1.,maximum_temperature_ev=20000.)]))
    fixture_rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,
        kind='power_law',A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=.7,
        minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    fixture=RecombiningPlasma(fixture_ion,fixture_rr,ionization_enabled=False)
    for steps in (50,100,200):
        result=evolve(fixture,1e20,[.2,.8],2.,.03,steps)
        x,t,p=power_law_reference(result['time_s'])
        result.update(reference_charge=x,reference_temperature_ev=t,reference_photons_ev_per_nucleus=p)
        np.savez_compressed(out/f'analytic_rr_t{steps}.npz',**result)
        entry=dict(steps=steps,maximum_fraction_error=float(abs(result['charge']-x).max()),
            maximum_temperature_error_ev=float(abs(result['temperature_ev']-t).max()),
            maximum_photon_energy_error_ev_per_nucleus=float(abs(result['state'][:,-1]/EV-p).max()),
            maximum_energy_balance_ev_per_nucleus=float(abs(result['energy_balance_j_per_nucleus']).max()/EV))
        report['reference'].append(entry);print('analytic RR',entry,flush=True)
        axes[1,0].semilogy(result['time_s'][1:]*1000,abs(result['charge'][1:]-x[1:]),label=f'{steps} steps')
    density=8.961429541213399e19
    cases=[('singly_ionized',1.,0.,10.,3.5e-6,(400,800,1600)),
           ('partial_seed',.1,0.,10.,3.5e-6,(400,800,1600)),
           ('heated_seed',.1,20.,10.,3.5e-6,(400,800,1600)),
           ('cold_recombining',1.,0.,1.,1e-3,(800,1600,3200))]
    for case,seed,deposit,te,duration,levels in cases:
        f=np.zeros(19);f[0]=1-seed;f[1]=seed;previous=None
        control=ionization_only(argon,density,f,te,duration,levels[-1],deposit)
        np.savez_compressed(out/f'argon_{case}_ionization_only.npz',**control)
        for steps in levels:
            result=evolve(model,density,f,te,duration,steps,deposit)
            filename=f'argon_{case}_t{steps}.npz';np.savez_compressed(out/filename,**result)
            entry=dict(case=case,steps=steps,filename=filename,nucleus_density_m3=density,
                initial_singly_ionized_fraction=seed,initial_temperature_ev=te,heating_ev_per_nucleus=deposit,
                duration_s=duration,heating_duration_s=1e-6,
                final_mean_charge=float(result['charge'][-1]),final_temperature_ev=float(result['temperature_ev'][-1]),
                final_photons_ev_per_nucleus=float(result['state'][-1,-1]/EV),
                final_populations=result['state'][-1,:-2].tolist(),minimum_population=float(result['state'][:,:-2].min()),
                maximum_number_error=float(abs(result['state'][:,:-2].sum(axis=1)-1).max()),
                maximum_energy_balance_ev_per_nucleus=float(abs(result['energy_balance_j_per_nucleus']).max()/EV),
                maximum_unknown_rr_range_event_fraction=float(result['unknown_rr_range_event_fraction'].max()),
                maximum_outside_known_rr_range_event_fraction=float(result['outside_known_rr_range_event_fraction'].max()),
                maximum_outside_ionization_fit_event_fraction=float(result['outside_ionization_fit_event_fraction'].max()),
                final_mean_charge_without_rr=float(control['charge'][-1]),
                final_temperature_without_rr_ev=float(control['temperature_ev'][-1]))
            if previous is not None:
                entry['maximum_fraction_change_from_double_dt']=float(abs(previous['state'][:,:-2]-result['state'][::2,:-2]).max())
                entry['maximum_temperature_change_from_double_dt_ev']=float(abs(previous['temperature_ev']-result['temperature_ev'][::2]).max())
                entry['maximum_photon_change_from_double_dt_ev']=float(abs(previous['state'][:,-1]-result['state'][::2,-1]).max()/EV)
            report['argon'].append(entry);previous=result;save()
            print('argon',case,steps,entry['final_mean_charge'],entry['final_temperature_ev'],entry['final_photons_ev_per_nucleus'],flush=True)
        label=case.replace('_',' ')
        if case=='cold_recombining':
            axes[1,1].plot(result['time_s']*1000,result['charge'],label='Mean charge')
            axes[1,1].plot(result['time_s']*1000,result['temperature_ev']/te,label='Te / initial Te')
        else:
            t=result['time_s']*1e6
            axes[0,0].plot(t,result['temperature_ev'],label=label)
            axes[0,1].plot(t,result['charge'],label=label)
            axes[0,2].plot(t,result['state'][:,-1]/EV*1000,label=label)
    labels=[('Argon electron temperature','Time [µs]','Te [eV]'),('Argon evolving charge','Time [µs]','Electrons per nucleus'),
        ('Escaping RR photon energy','Time [µs]','Photon energy [meV/nucleus]'),
        ('Independent recombination solution','Time [ms]','Absolute ion-fraction error'),
        ('1 eV source-only recombination','Time [ms]','Relative temperature / mean charge'),
        ('Published original-rate comparison','Temperature [eV]','Fit / original rate − 1 [%]')]
    for ax,(title,xlabel,ylabel) in zip(axes.flat,labels):
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2 if ax==axes[1,2] else 1)
    fig.suptitle('Ionization and radiative recombination evolve charge, electron energy and escaping light\nFixed-volume component studies; instantaneous ground-state cascades, no opacity or other atomic channels')
    fig.savefig(out/'recombination_evolution.png',dpi=170);fig.savefig(out/'recombination_evolution.pdf');plt.close(fig)
    report['terminal']=True;save();print('Recombination source calculation complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination')
    run(parser.parse_args().destination)
