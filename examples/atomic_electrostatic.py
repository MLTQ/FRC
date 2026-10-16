"""Spatial two-fluid plasma waves with charge-continuous electric work."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV, epsilon_0 as EPS0
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_electrostatic import AtomicElectrostatic
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh spatial electrostatic output required')
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    n=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=False,recombination_enabled=False,radiation=False)
    source=FixedFrameAtomicMoments(n,39.948*atomic_mass,[0.,0.,0.])
    density=1e14;length=.1;gamma=5/3;temperature=np.array([.2,2.]);masses=source.mass[[1,-1]]
    k=2*np.pi/length;plasma2=EV**2*density/(EPS0*masses);sound2=gamma*EV*temperature/masses
    a,b=k*k*sound2+plasma2
    frequency2=.5*(a+b+np.sqrt((a-b)**2+4*plasma2.prod()))
    frequency=np.sqrt(frequency2);ratio=plasma2[0]/(a-frequency2);duration=2.3*2*np.pi/frequency
    identity=provenance(['examples.atomic_electrostatic']);out.mkdir(parents=True);freeze(out,identity)
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],
        reference_density_m3=density,length_m=length,temperature_ev=temperature.tolist(),frequency_rad_s=float(frequency),
        ion_to_electron_density_amplitude=float(ratio),duration_s=float(duration),
        scope='Complete nonlinear Euler-Poisson standing-wave histories against an independent linear two-fluid dispersion relation. Reflecting verification interval, physical masses, no reactions/magnetic/free-boundary or formation claim.')
    began=time.monotonic();histories={}
    def save():
        report['elapsed_seconds']=time.monotonic()-began;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    for amplitude in (1e-3,1e-4):
        amplitudes=amplitude*np.array([ratio,1.])
        for cells in (16,32,64):
            model=AtomicElectrostatic(source,cells,length,density,maximum_plasma_phase=.04*16/cells)
            profile=np.cos(k*model.x)*np.sinc(1/cells);face_x=np.arange(cells+1)*model.dx
            numbers=1+profile[:,None]*amplitudes;pressure=temperature*(1+gamma*profile[:,None]*amplitudes)
            initial=[]
            for f,p in zip(numbers,pressure):
                atom=n.initial([0.,f[0]],p[1]/f[1]);atom[n.electrons]=f[1];atom[n.thermal]=1.5*p[1]
                initial.append(source.initial(atom,p[0]/f[0],np.zeros((3,3))))
            state=np.asarray(initial);t=0.;times=[0.];states=[state.copy()];failure=None
            max_ampere=0.;max_work=0.
            while t<duration:
                dt=min(model.timestep(state),duration-t)
                try:candidate=model.step(state,t,dt,reactions=False)
                except (ValueError,RuntimeError) as error:
                    failure=dict(step=len(times),time_s=t,dt_s=dt,kind=type(error).__name__,message=str(error));break
                state=candidate;t=min(t+dt,duration);times.append(t);states.append(state.copy())
                max_ampere=max(max_ampere,model.last_diagnostic['ampere_residual_v_m_s'])
                max_work=max(max_work,model.last_diagnostic['electric_work_residual_ev_m_s'])
            times=np.asarray(times);states=np.asarray(states);fields=np.asarray([model.electric(s) for s in states])
            rows=[model.totals(s) for s in states];budgets={key:np.asarray([r[key] for r in rows]) for key in rows[0]}
            exact_number=1+np.cos(frequency*times)[:,None,None]*profile[None,:,None]*amplitudes
            field_amplitude=EV*density*(amplitudes[0]-amplitudes[1])/(EPS0*k)
            exact_electric=field_amplitude*np.cos(frequency*times)[:,None]*np.sin(k*face_x)
            exact_current=EV*density*(amplitudes[0]-amplitudes[1])*frequency/k*np.sin(frequency*times)[:,None]*np.sin(k*model.x)*np.sinc(1/cells)
            currents=EV*density*source.momentum_scale*np.einsum('s,tcs->tc',np.array([0.,1.,-1.])/source.mass,states[:,:,source.momentum[:,2]])
            particle_error=abs(states[:,:,[n.ground[1],n.electrons]]-exact_number)
            case=dict(name=f'wave_a{amplitude:g}_n{cells}',amplitude=amplitude,cells=cells,steps=len(times)-1,
                success=failure is None,failure=failure,completed_time_s=float(times[-1]),maximum_dt_s=float(np.diff(times).max(initial=0.)),
                maximum_electron_density_error_over_amplitude=float(particle_error[:,:,1].max()/amplitude),
                final_electron_density_mean_error_over_amplitude=float(particle_error[-1,:,1].mean()/amplitude),
                maximum_electric_error_over_initial_amplitude=float(abs(fields-exact_electric).max()/abs(field_amplitude)),
                maximum_current_error_over_mode_amplitude=float(abs(currents-exact_current).max()/(abs(field_amplitude)*EPS0*frequency)),
                maximum_nucleus_integral_error=float(abs(budgets['nuclei']-budgets['nuclei'][0]).max()),
                maximum_charge_integral_error=float(abs(budgets['charge']-budgets['charge'][0]).max()),
                maximum_total_energy_error_ev_m=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()),
                maximum_linear_field_energy_error_ev_m=float(abs(budgets['linear_field_total_energy_ev']-budgets['linear_field_total_energy_ev'][0]).max()),
                maximum_total_energy_error_over_initial_field_energy=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()/budgets['electric_energy_ev'][0]),
                maximum_gauss_residual_over_reference_charge=float(budgets['maximum_gauss_residual_over_reference_charge'].max()),
                maximum_stage_ampere_residual_over_mode_rate=float(max_ampere/(abs(field_amplitude)*frequency)),
                maximum_stage_electric_work_residual_ev_m_s=float(max_work),
                maximum_neutral_number=float(states[:,:,0].max()))
            history=dict(time_s=times,states=states,x_m=model.x,face_x_m=face_x,electric_v_m=fields,
                analytic_electric_v_m=exact_electric,analytic_number=exact_number,current_a_m2=currents,analytic_current_a_m2=exact_current,**budgets)
            np.savez_compressed(out/(case['name']+'.npz'),**history);histories[amplitude,cells]=history
            report['cases'].append(case);save();print('Spatial electrostatic wave',case,flush=True)
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for cells in (16,32,64):
        h=histories[1e-4,cells];axes[0,0].plot(h['x_m'],(h['states'][-1,:,n.electrons]-1)/1e-4,label=f'{cells} cells')
        axes[0,1].plot(h['face_x_m'],h['electric_v_m'][-1],label=f'{cells} cells')
    axes[0,0].plot(h['x_m'],(h['analytic_number'][-1,:,1]-1)/1e-4,'k--',label='Linear two-fluid mode')
    axes[0,1].plot(h['face_x_m'],h['analytic_electric_v_m'][-1],'k--',label='Linear Gauss field')
    for amplitude in (1e-3,1e-4):
        rows=[r for r in report['cases'] if r['amplitude']==amplitude]
        axes[1,0].loglog([r['cells'] for r in rows],[r['maximum_electric_error_over_initial_amplitude'] for r in rows],'o-',label=f'Amplitude {amplitude:g}')
        h=histories[amplitude,64]
        axes[1,1].plot(h['time_s']*1e9,(h['total_energy_ev']-h['total_energy_ev'][0])/h['electric_energy_ev'][0],label=f'Amplitude {amplitude:g}')
    for ax,xlabel,ylabel in zip(axes.flat,('Axial coordinate [m]','Axial coordinate [m]','Cells','Time [ns]'),('Electron perturbation / amplitude','Final electric field [V/m]','Maximum E error / initial amplitude','Total energy drift / initial field energy')):
        ax.set(xlabel=xlabel,ylabel=ylabel);ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Spatial species/field feedback: full Euler-Poisson histories\nReflecting interval; analytic linear-mode comparison, no formation claim')
    fig.savefig(out/'spatial_electrostatic.png',dpi=170);fig.savefig(out/'spatial_electrostatic.pdf');plt.close(fig)
    if provenance(['examples.atomic_electrostatic'])!=identity:raise RuntimeError('Source changed during spatial field evolution')
    report['terminal']=True;report['success']=all(r['success'] for r in report['cases']);save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');a=p.parse_args()
    r=run(a.destination);raise SystemExit(0 if r['success'] else 1)
