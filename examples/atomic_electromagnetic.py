"""Complete reacting field/current histories and independent uniform modes."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV, epsilon_0 as EPS0
from plasma_bench.atomic_electromagnetic import AtomicElectromagnetic
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_collisions import CollisionalAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .atomic_network import argon_models
from .implicit_plasma_records import provenance, freeze


def run(destination,neutral_path,ion_path,ionization_path,*,method='midpoint'):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh electromagnetic-source output required')
    models,inputs=argon_models(neutral_path,ion_path,ionization_path)
    identity=provenance(['examples.atomic_electromagnetic']);out.mkdir(parents=True);freeze(out,identity)
    report=dict(terminal=False,source_sha256=identity['source_sha256'],inputs=inputs,cases=[],electromagnetic_method=method,
        scope='Homogeneous multifluid-Maxwell source with full reacting/colliding species. Physical masses, fixed laboratory B, no spatial curls/boundary or formation claim. Midpoint stability does not resolve plasma/cyclotron phase.')
    started=time.monotonic();histories={}
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def record(name,model,initial,electric,duration,steps,*,sources,heat=None,reference=None):
        times=np.linspace(0.,duration,steps+1);state=initial.copy();field=np.asarray(electric).copy()
        states=[state.copy()];fields=[field.copy()];failure=None
        for k,(a,b) in enumerate(zip(times[:-1],times[1:]),1):
            try:candidate,next_field=model.step(state,field,a,b-a,sources=sources,heat=heat)
            except (ValueError,RuntimeError) as error:
                failure=dict(step=k,time_s=a,dt_s=b-a,kind=type(error).__name__,message=str(error));break
            state,field=candidate,next_field;states.append(state.copy());fields.append(field.copy())
        states=np.asarray(states);fields=np.asarray(fields);times=times[:len(states)]
        rows=[model.observe(s,e) for s,e in zip(states,fields)]
        history={key:np.asarray([row[key] for row in rows]) for key in rows[0]}
        history.update(time_s=times,states=states,reference_velocity_m_s=model.frame,magnetic_field_t=model.magnetic)
        case=dict(name=name,success=failure is None,failure=failure,requested_steps=steps,accepted_steps=len(states)-1,
            duration_s=duration,completed_time_s=float(times[-1]),dt_s=duration/steps,
            reference_density_m3=model.density,magnetic_field_t=model.magnetic.tolist(),
            maximum_plasma_frequency_times_dt=float(history['plasma_frequency_rad_s'].max()*duration/steps),
            maximum_cyclotron_frequency_times_dt=float(history['maximum_cyclotron_frequency_rad_s'].max()*duration/steps),
            maximum_energy_budget_error_ev=float(abs(history['total_energy_budget_ev']-history['total_energy_budget_ev'][0]).max()),
            maximum_total_momentum_error_scaled=float(np.linalg.norm(history['total_momentum_scaled']-history['total_momentum_scaled'][0],axis=1).max()),
            maximum_charge_error=float(abs(history['charge_balance']-history['charge_balance'][0]).max()),
            maximum_nucleus_error=float(abs(history['number']-history['number'][0]).max()),
            final_electrons=float(history['electrons'][-1]),final_electron_temperature_ev=float(history['electron_temperature_ev'][-1]),
            final_current_a_m2=history['current_a_m2'][-1].tolist(),final_electric_v_m=fields[-1].tolist())
        if reference is not None:
            exact_e,exact_v=reference(times)
            f=[model.source.fields(s) for s in states];velocity=np.asarray([x['velocity']+model.frame for x in f])
            exact_j=model.density*EV*(exact_v[:,1]-exact_v[:,-1])
            history.update(analytic_electric_v_m=exact_e,analytic_velocity_m_s=exact_v,velocity_m_s=velocity)
            case.update(maximum_electric_error_v_m=float(np.linalg.norm(fields-exact_e,axis=1).max()),
                maximum_current_error_a_m2=float(np.linalg.norm(history['current_a_m2']-exact_j,axis=1).max()),
                maximum_occupied_species_velocity_error_m_s=float(np.linalg.norm(velocity[:,1:]-exact_v[:,1:],axis=2).max()),
                maximum_species_temperature_error_ev=float(abs(np.asarray([x['temperature_ev'] for x in f])-f[0]['temperature_ev']).max()))
        np.savez_compressed(out/f'{name}.npz',**history)
        histories[name]=history;report['cases'].append(case);save();print('Atomic EM',case,flush=True)
    save()
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    n=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=False,recombination_enabled=False,radiation=False)
    source=FixedFrameAtomicMoments(n,39.948*atomic_mass,[200.,-100.,300.])
    density=1e16;field=100.;magnetic=.02;charges=np.array([EV,-EV]);masses=source.mass[[1,-1]]
    plasma2=density*charges**2/(EPS0*masses);frequency=np.sqrt(plasma2.sum())
    initial=source.initial(n.initial([0.,1.],2.),.2,np.zeros((n.count+1,3)))
    model=AtomicElectromagnetic(source,density,[0.,0.,magnetic],method=method)
    def parallel(t):
        electric=np.zeros((len(t),3));electric[:,2]=field*np.cos(frequency*t)
        velocity=np.zeros((len(t),3,3));velocity[:,1:,2]=np.sin(frequency*t)[:,None]*(charges/masses*field/frequency)
        return electric,velocity
    # Uniform perpendicular cold-fluid modes: omega=0 and the two nonzero
    # roots of omega=sum(omega_ps^2/(omega-Omega_s)). No numerical integrator.
    cyclotron=charges*magnetic/masses
    spread=np.sqrt((cyclotron[0]-cyclotron[1])**2+4*plasma2.sum())
    frequencies=np.r_[0.,.5*(cyclotron.sum()+spread),.5*(cyclotron.sum()-spread)]
    eigenvectors=np.vstack((np.ones(3),1j*(charges/masses)[:,None]/(frequencies[None,:]-cyclotron[:,None])))
    coefficients=np.linalg.solve(eigenvectors,np.array([field,0.,0.]))
    def perpendicular(t):
        modes=(np.exp(-1j*t[:,None]*frequencies)*coefficients)@eigenvectors.T
        electric=np.zeros((len(t),3));electric[:,:2]=np.column_stack((modes[:,0].real,modes[:,0].imag))
        velocity=np.zeros((len(t),3,3));velocity[:,1:,:2]=np.stack((modes[:,1:].real,modes[:,1:].imag),axis=-1)
        return electric,velocity
    for name,reference,duration,electric in [
            ('parallel',parallel,5.3*2*np.pi/frequency,[0.,0.,field]),
            ('perpendicular',perpendicular,5.3*2*np.pi/max(abs(frequencies)),[field,0.,0.])]:
        for steps in (100,200,400):record(f'{name}_{steps}',model,initial,electric,duration,steps,sources=False,reference=reference)
    n=models['joint'];source=CollisionalAtomicMoments(n,39.948*atomic_mass,[0.,0.,0.],coulomb_log=10.)
    velocity=np.zeros((n.count+1,3));velocity[:-1,2]=np.where(n.charge==0,1000.,-1000.);velocity[-1]=[0.,10000.,0.]
    initial=source.initial(n.initial(np.r_[.9,.1,np.zeros(17)],10.),.1,velocity)
    model=AtomicElectromagnetic(source,8.961429541213399e19,[0.,0.,.02],method=method)
    def heat(t):return 40./1e-6*np.sin(np.pi*t/1e-6)**2 if 0<t<1e-6 else 0.
    for steps in (700,1400):record(f'heated_argon_{steps}',model,initial,[0.,0.,0.],3.5e-6,steps,sources=True,heat=heat)
    coarse,fine=histories['heated_argon_700'],histories['heated_argon_1400']
    if len(coarse['states'])==701 and len(fine['states'])==1401:
        report['reacting_time_halving']=dict(
            maximum_population_l1=float(abs(coarse['states'][:,:n.count]-fine['states'][::2,:n.count]).sum(axis=1).max()),
            maximum_electric_difference_v_m=float(np.linalg.norm(coarse['electric_v_m']-fine['electric_v_m'][::2],axis=1).max()),
            maximum_current_difference_a_m2=float(np.linalg.norm(coarse['current_a_m2']-fine['current_a_m2'][::2],axis=1).max()),
            maximum_electron_temperature_difference_ev=float(abs(coarse['electron_temperature_ev']-fine['electron_temperature_ev'][::2]).max()),
            scope='Full coupled-source time sensitivity. Large omega_p*dt does not resolve fast phase even if integral budgets are small.')
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    h=histories['perpendicular_400'];t=h['time_s']*1e9
    axes[0,0].plot(t,h['electric_v_m'][:,0],label='Computed Ex');axes[0,0].plot(t,h['analytic_electric_v_m'][:,0],'k--',label='Analytic cold-fluid modes')
    for ax,key in [(axes[0,1],'electrons'),(axes[1,0],'current_a_m2'),(axes[1,1],'total_energy_budget_ev')]:
        for steps in (700,1400):
            h=histories[f'heated_argon_{steps}'];y=h[key]
            if key=='current_a_m2':y=np.linalg.norm(y,axis=1)
            if key=='total_energy_budget_ev':y=y-y[0]
            ax.plot(h['time_s']*1e6,y,label=f'{steps} requested steps'+(' (incomplete)' if len(h['states'])!=steps+1 else ''))
    for ax,label,xlabel in zip(axes.flat,('Perpendicular Ex [V/m]','Electrons / nucleus','Current magnitude [A/m²]','Particle + field energy budget [eV/nucleus]'),('Time [ns]','Time [µs]','Time [µs]','Time [µs]')):
        ax.set(xlabel=xlabel,ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(f'Reacting species with electric-field feedback ({method})\nHomogeneous source; phase accuracy remains separate from energy conservation')
    fig.savefig(out/'reacting_electromagnetic.png',dpi=170);fig.savefig(out/'reacting_electromagnetic.pdf');plt.close(fig)
    if provenance(['examples.atomic_electromagnetic'])!=identity:raise RuntimeError('Source changed during electromagnetic history')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ('destination','neutral_input','ion_input','ionization_input'):p.add_argument(key)
    p.add_argument('--method',choices=('midpoint','exponential'),default='midpoint')
    a=p.parse_args();r=run(a.destination,a.neutral_input,a.ion_input,a.ionization_input,method=a.method)
    raise SystemExit(0 if r['success'] else 1)
