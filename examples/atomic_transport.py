"""Full spatial species histories against translating contact and decay laws."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_transport import AtomicTransport
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .atomic_network import argon_models
from .implicit_plasma_records import provenance, freeze


def run(destination, neutral_path, ion_path, ionization_path):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh spatial-species output directory required')
    models,inputs=argon_models(neutral_path,ion_path,ionization_path)
    identity=provenance(['examples.atomic_transport']);out.mkdir(parents=True);freeze(out,identity)
    started=time.monotonic();report=dict(terminal=False,source_sha256=identity['source_sha256'],inputs=inputs,cases=[],
        scope='Periodic axial species-transport components. Exact zero-current plasma contact and exactly uncharged translating neutral decay; no arbitrary charged-field, lab-end or formation claim.')
    finals={}
    def save():
        report['elapsed_seconds']=time.monotonic()-started;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    for kind,cells_list in [('contact',(32,64,128)),('neutral_decay',(16,32,64))]:
        if kind=='contact':
            network=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=False,recombination_enabled=False,radiation=False)
            length,speed,temperature,duration=1.,1e5,.2,1e-5
        else:
            joint=models['joint']
            network=AtomicNetwork(joint.ionization,joint.recombination,{0:joint.levels[0]},
                ionization_enabled=False,recombination_enabled=False)
            length,speed,temperature,duration=.01,2.5e5,.02585,40e-9
            radiation=network.levels[0].radiation;upper=int(radiation.upper[0]);lower=int(radiation.lower[0])
            if np.count_nonzero(radiation.upper==upper)!=1 or np.any(radiation.upper==lower):
                raise ValueError('Independent decay requires one radiative destination and stable lower level')
            decay=float(radiation.A[0]);energy=network.energy_ev[upper]-network.energy_ev[lower]
        model=FixedFrameAtomicMoments(network,39.948*atomic_mass,[0.,0.,speed])
        velocity=np.tile([0.,0.,speed],(network.count+1,1))
        for cells in cells_list:
            transport=AtomicTransport(model,cells,length,1e19)
            amplitude=.2*np.sinc(1/cells);phase=2*np.pi*transport.x/length;omega=2*np.pi*speed/length
            profile=1+amplitude*np.sin(phase);initial=[]
            for number in profile:
                fractions=np.zeros(len(network.ground));fractions[1 if kind=='contact' else 0]=number
                atomic=network.initial(fractions,temperature/number if kind=='contact' else 0.)
                if kind=='neutral_decay':atomic[network.ground[0]]=0.;atomic[upper]=number
                initial.append(model.initial(atomic,temperature/number,velocity))
            initial=np.asarray(initial)
            def exact(t):
                density=1+amplitude*np.sin(phase-omega*t)
                result=initial.copy();result[:,:]=0.
                if kind=='contact':
                    result[:,network.ground[1]]=density;result[:,network.electrons]=density
                    result[:,network.thermal]=1.5*temperature;result[:,model.heavy_energy[network.ground[1]]]=1.5*temperature
                else:
                    remaining=np.exp(-decay*t)
                    result[:,upper]=remaining*density;result[:,lower]=(1-remaining)*density
                    result[:,model.heavy_energy[upper]]=remaining*1.5*temperature
                    result[:,model.heavy_energy[lower]]=(1-remaining)*1.5*temperature
                    integral=(1-remaining)+amplitude*np.imag(np.exp(1j*phase)*decay*(-np.expm1(-(decay+1j*omega)*t))/(decay+1j*omega))
                    result[:,network.photons]=energy*integral
                return result
            state=initial.copy();t=0.;times=[0.];states=[state.copy()];failure=None
            while t<duration:
                dt=min(transport.timestep(state),duration-t)
                try:candidate=transport.step(state,t,dt,reactions=kind=='neutral_decay')
                except (ValueError,RuntimeError) as error:
                    failure=dict(step=len(times),time_s=t,dt_s=dt,kind=type(error).__name__,message=str(error));break
                state=candidate;t=min(t+dt,duration);times.append(t);states.append(state.copy())
            times=np.asarray(times);states=np.asarray(states);rows=[transport.totals(s) for s in states]
            budgets={key:np.asarray([row[key] for row in rows]) for key in rows[0]}
            reference=exact(times[-1]);delta=state-reference
            initial_number=budgets['nuclei'][0]
            charge=states[:,:,network.electrons]-np.einsum('s,tcs->tc',network.charge,states[:,:,:network.count])
            case=dict(name=f'{kind}_{cells}',kind=kind,cells=cells,steps=len(times)-1,success=failure is None,failure=failure,
                duration_s=duration,completed_time_s=float(times[-1]),dx_m=transport.dx,maximum_dt_s=float(np.diff(times).max(initial=0.)),
                final_population_mean_l1=float(abs(delta[:,:network.count]).sum(axis=1).mean()),
                final_population_max_error=float(abs(delta[:,:network.count]).max()),
                final_photon_max_error_ev=float(abs(delta[:,network.photons]).max()),
                maximum_nucleus_integral_error=float(abs(budgets['nuclei']-initial_number).max()),
                maximum_mass_integral_error=float(abs(budgets['mass_over_neutral']-budgets['mass_over_neutral'][0]).max()),
                maximum_charge_integral_error=float(abs(budgets['charge']).max()),maximum_local_charge_error=float(abs(charge).max()),
                maximum_momentum_integral_error=float(np.linalg.norm(budgets['momentum_scaled']-budgets['momentum_scaled'][0],axis=1).max()),
                maximum_energy_integral_error_ev_m=float(abs(budgets['energy_ev']-budgets['energy_ev'][0]).max()),
                maximum_electron_number=float(states[:,:,network.electrons].max()))
            np.savez_compressed(out/f'{kind}_{cells}.npz',time_s=times,x_m=transport.x,states=states,analytic_final=reference,
                reference_velocity_m_s=model.reference_velocity,**budgets)
            report['cases'].append(case);finals[kind,cells]=(transport.x,state,reference,network);save()
            print('Atomic transport',case,flush=True)
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for kind,cells_list,ax in [('contact',(32,64,128),axes[0,0]),('neutral_decay',(16,32,64),axes[0,1])]:
        for cells in cells_list:
            x,state,reference,n=finals[kind,cells]
            ax.plot(x,state[:,:n.count].sum(axis=1),label=f'{cells} cells')
        ax.plot(x,reference[:,:n.count].sum(axis=1),'k--',label='Exact cell average')
        ax.set(title=kind.replace('_',' '),xlabel='Axial coordinate [m]',ylabel='Nuclei / reference density');ax.legend();ax.grid(alpha=.2)
    for kind in ('contact','neutral_decay'):
        cases=[c for c in report['cases'] if c['kind']==kind and c['success']]
        axes[1,0].loglog([c['cells'] for c in cases],[c['final_population_mean_l1'] for c in cases],'o-',label=kind)
    axes[1,0].set(xlabel='Cells',ylabel='Final mean population L1 error');axes[1,0].legend();axes[1,0].grid(alpha=.2)
    x,state,reference,n=finals['neutral_decay',64]
    axes[1,1].plot(x,state[:,n.photons],label='Produced locally');axes[1,1].plot(x,reference[:,n.photons],'k--',label='Independent decay/advection integral')
    axes[1,1].set(xlabel='Axial coordinate [m]',ylabel='Produced photons [eV / reference nucleus]');axes[1,1].legend();axes[1,1].grid(alpha=.2)
    fig.suptitle('Spatial species transport and local decay — periodic mathematical limits; no field feedback')
    fig.savefig(out/'spatial_species.png',dpi=170);fig.savefig(out/'spatial_species.pdf');plt.close(fig)
    if provenance(['examples.atomic_transport'])!=identity:raise RuntimeError('Source changed during spatial evolution')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save()
    print('Atomic transport terminal',report['success'],out,flush=True);return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('destination','neutral_input','ion_input','ionization_input'):parser.add_argument(name)
    a=parser.parse_args();report=run(a.destination,a.neutral_input,a.ion_input,a.ionization_input)
    raise SystemExit(0 if report['success'] else 1)
