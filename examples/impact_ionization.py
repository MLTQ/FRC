"""Complete direct-ionization source evolutions with independent scalar reference."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV
from scipy.integrate import quad
from scipy.optimize import brentq
from plasma_bench.impact_ionization import ImpactIonization
from .implicit_plasma_records import provenance, freeze


def evolve(model, density, fraction, temperature_ev, duration, steps, deposit_ev=0., pulse=1e-6):
    f = np.asarray(fraction, float)
    y = np.r_[f, 1.5*(model.charge@f)*temperature_ev*EV]
    states = [y.copy()]; times = np.linspace(0., duration, steps+1)
    def rhs(t, state):
        heat = 2*deposit_ev*EV/pulse*np.sin(np.pi*t/pulse)**2 if 0<t<pulse else 0.
        df, de, _ = model.source(state[:-1], state[-1], density, heating_j_per_nucleus_s=heat)
        return np.r_[df, de]
    for t, end in zip(times[:-1], times[1:]):
        h = end-t
        k1 = rhs(t,y); k2 = rhs(t+h/2,y+h*k1/2)
        k3 = rhs(t+h/2,y+h*k2/2); k4 = rhs(end,y+h*k3)
        y = y+h*(k1+2*k2+2*k3+k4)/6
        model.source(y[:-1],y[-1],density)  # Domain check, never a state repair.
        states.append(y.copy())
    states = np.asarray(states); populations = states[:,:-1]
    charge = populations@model.charge
    temperature = states[:,-1]/(1.5*charge*EV)
    atomic = populations@model.atomic_energy
    inside = np.minimum(times/pulse,1.)
    heating = deposit_ev*EV*(inside-np.sin(2*np.pi*inside)/(2*np.pi))
    total = atomic+states[:,-1]
    diagnostic = [model.source(f,e,density)[2] for f,e in zip(populations,states[:,-1])]
    return dict(time_s=times,state=states,charge=charge,temperature_ev=temperature,
        atomic_energy_j_per_nucleus=atomic,analytic_heating_j_per_nucleus=heating,
        energy_balance_j_per_nucleus=total-total[0]-heating,
        outside_fit_event_fraction=np.array([d['outside_fit_event_fraction'] or 0. for d in diagnostic]),
        discrepant_threshold_event_fraction=np.array([d['discrepant_threshold_event_fraction'] or 0. for d in diagnostic]))


def hydrogen_reference(times, density, x0, t0):
    # Direct analytic quadrature, independent of production source/rate functions.
    chi=13.6; total=1.5*x0*t0+chi*x0
    def temperature(x):
        return (total-chi*x)/(1.5*x)
    def inverse_rate(x):
        u=13.6/temperature(x)
        k=2.91e-14*u**.39*np.exp(-u)/(.232+u)
        return 1/(density*k*x*(1-x))
    upper=min(.999999,total/(chi+1.5))  # Reference stays at Te>=1eV.
    def elapsed(x):
        return quad(inverse_rate,x0,x,epsabs=1e-16,epsrel=2e-12)[0]
    if elapsed(upper)<times[-1]:
        raise ValueError('Requested reference exceeds its Te>=1eV interval')
    ion=np.array([x0 if t==0 else brentq(lambda x:elapsed(x)-t,x0,upper,xtol=5e-15) for t in times])
    return ion,temperature(ion)


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('A new output directory is required')
    out.mkdir(parents=True);start=time.monotonic()
    identity=provenance(['examples.impact_ionization']);freeze(out,identity)
    data=Path(__file__).resolve().parents[1]/'data/atomic/argon'
    shutil.copytree(data,out/'atomic_inputs')
    argon=ImpactIonization.from_file(data/'argon.json')
    hydrogen=ImpactIonization(dict(nuclear_charge=1,transitions=[dict(charge=0,
        threshold_fit_ev=13.6,ionization_ev=13.6,A_cm3_s=2.91e-8,P=0.,X=.232,K=.39,
        minimum_temperature_ev=1.,maximum_temperature_ev=20000.)]))
    report=dict(terminal=False,source_sha256=identity['source_sha256'],atomic_input_sha256=argon.input_sha256,
        scope='Homogeneous direct-ionization source only; no recombination, excitation, radiation, metastables, motion, transport or preparation prediction.',
        hydrogen_reference='Independent t(x) quadrature with energy-derived Te(x), N=1e20m^-3, x0=.2, Te0=30eV, chi=13.6eV.',
        hydrogen=[],argon=[])
    def save():
        report['elapsed_seconds']=time.monotonic()-start
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for steps in (100,200,400):
        result=evolve(hydrogen,1e20,[.8,.2],30.,20e-6,steps)
        reference,temp=hydrogen_reference(result['time_s'],1e20,.2,30.)
        result['reference_ion_fraction']=reference;result['reference_temperature_ev']=temp
        filename=f'hydrogen_t{steps}.npz';np.savez_compressed(out/filename,**result)
        row=dict(steps=steps,filename=filename,maximum_fraction_error=float(abs(result['charge']-reference).max()),
            maximum_temperature_error_ev=float(abs(result['temperature_ev']-temp).max()),
            maximum_energy_balance_ev_per_nucleus=float(abs(result['energy_balance_j_per_nucleus']).max()/EV))
        report['hydrogen'].append(row);print('hydrogen',row,flush=True)
        axes[1,0].semilogy(result['time_s'][1:]*1e6,abs(result['charge'][1:]-reference[1:]),label=f'{steps} steps')
    density=8.961429541213399e19
    for case,seed,deposit in [('singly_ionized',1.,0.),('partial_seed',.1,0.),('heated_seed',.1,20.)]:
        f=np.zeros(19);f[0]=1-seed;f[1]=seed
        previous=None
        for steps in (400,800,1600):
            result=evolve(argon,density,f,10.,3.5e-6,steps,deposit)
            filename=f'argon_{case}_t{steps}.npz';np.savez_compressed(out/filename,**result)
            row=dict(case=case,steps=steps,filename=filename,nucleus_density_m3=density,
                initial_singly_ionized_fraction=seed,initial_temperature_ev=10.,heating_ev_per_nucleus=deposit,
                duration_s=3.5e-6,heating_duration_s=1e-6,
                final_mean_charge=float(result['charge'][-1]),final_temperature_ev=float(result['temperature_ev'][-1]),
                final_populations=result['state'][-1,:-1].tolist(),minimum_population=float(result['state'][:,:-1].min()),
                maximum_number_error=float(abs(result['state'][:,:-1].sum(axis=1)-1).max()),
                maximum_energy_balance_ev_per_nucleus=float(abs(result['energy_balance_j_per_nucleus']).max()/EV),
                maximum_outside_fit_event_fraction=float(result['outside_fit_event_fraction'].max()),
                maximum_discrepant_threshold_event_fraction=float(result['discrepant_threshold_event_fraction'].max()))
            if previous is not None:
                row['maximum_fraction_change_from_double_dt']=float(abs(previous['state'][:,:-1]-result['state'][::2,:-1]).max())
                row['maximum_temperature_change_from_double_dt_ev']=float(abs(previous['temperature_ev']-result['temperature_ev'][::2]).max())
            report['argon'].append(row);previous=result
            print('argon',case,steps,row['final_mean_charge'],row['final_temperature_ev'],flush=True)
        t=result['time_s']*1e6
        axes[0,0].plot(t,result['temperature_ev'],label=case)
        axes[0,1].plot(t,result['charge'],label=case)
        axes[0,2].plot(t,(result['atomic_energy_j_per_nucleus']-result['atomic_energy_j_per_nucleus'][0])/EV,label=case)
        axes[1,1].plot(t,result['energy_balance_j_per_nucleus']/EV,label=case)
        axes[1,2].plot(t,result['outside_fit_event_fraction'],label=case)
        save()
    labels=[('Argon electron temperature','Te [eV]'),('Argon evolving charge','Electrons per nucleus'),
        ('Energy spent on further ionization','Atomic energy gain [eV/nucleus]'),
        ('Independent hydrogen quadrature','Absolute ion-fraction error'),
        ('Thermal + atomic energy − heating','Energy residual [eV/nucleus]'),
        ('Rate-fit extrapolation content','Fraction of ionization events')]
    for ax,(title,ylabel) in zip(axes.flat,labels):
        ax.set(title=title,xlabel='Time [µs]',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Electron-impact ionization consumes thermal energy and changes electron density\nSource-only fixed-volume studies; ground-state rates, no recombination/radiation/transport')
    fig.savefig(out/'ionization_evolution.png',dpi=170);fig.savefig(out/'ionization_evolution.pdf');plt.close(fig)
    report['terminal']=True;save();print('Atomic source calculation complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination')
    run(parser.parse_args().destination)
