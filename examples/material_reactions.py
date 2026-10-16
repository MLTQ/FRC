"""Complete support-aware ionization source histories with analytic moments."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.species_material import MaterialSpecies
from plasma_bench.material_reactions import MaterialReactionOverlap,ReactionOverlapState
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance,freeze


def exact_moments(source,initial,time_s,density):
    """Independent local logistic number law and donor moment kinematics."""
    n=source.network;initial=np.asarray(initial);out=initial.copy()
    neutral=initial[...,0];electron=initial[...,n.electrons];reservoir=neutral+electron
    active=(neutral>0)&(electron>0);converted=np.zeros_like(neutral)
    updated=reservoir[active]/(1+(reservoir[active]/electron[active]-1)*np.exp(-1e-14*density*reservoir[active]*time_s))
    converted[active]=updated-electron[active]
    fraction=np.divide(converted,neutral,out=np.zeros_like(converted),where=neutral>0)
    out[...,0]-=converted;out[...,1]+=converted;out[...,n.electrons]+=converted
    parent_p=initial[...,source.momentum[0]];parent_e=initial[...,source.energy[0]]
    out[...,source.momentum[0]]-=fraction[...,None]*parent_p
    out[...,source.energy[0]]-=fraction*parent_e
    for child in (1,2):
        ratio=source.mass[child]/source.mass[0]
        out[...,source.momentum[child]]+=fraction[...,None]*ratio*parent_p
        out[...,source.energy[child]]+=fraction*ratio*parent_e
    out[...,n.thermal]-=13.6*converted
    return out


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh reaction-overlap output required')
    identity=provenance(['examples.material_reactions']);out.mkdir(parents=True);freeze(out,identity)
    began=time.monotonic();density=1e20;duration=2e-6;report=dict(terminal=False,cases=[],source_sha256=identity['source_sha256'],
        scope='Frozen-geometry atomic source substep from distinct physical material supports. Constant-coefficient single-ion mathematical limit with hydrogen masses/13.6eV binding; not a physical rate fit, moving-plasma evolution or formation prediction.')
    def save():
        report['elapsed_seconds']=time.monotonic()-began;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();histories={}
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',A_cm3_s=1e-12,
        reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    network=AtomicNetwork(ConstantIonization(),rr,recombination_enabled=False,radiation=False)
    boost=np.array([20000.,-10000.,5000.]);reference=FixedFrameAtomicMoments(network,1.00794*atomic_mass,boost)
    faces=np.array([np.linspace(-a,a,9) for a in (.015,.005,.010)])
    number=density*np.array([.8,.2,.2])[:,None]*np.diff(faces,axis=1)
    material=MaterialSpecies(number,reference.mass,[0.,1.,-1.],(-.03,.03))
    velocity=np.broadcast_to(boost,(3,9,3)).copy()+np.array([[200.,300.,100.],[-200.,400.,0.],[2000.,-1500.,-100.]])[:,None,:]
    velocity+=faces[...,None]/.015*np.array([1000.,200.,500.])
    initial_material=material.state(faces,np.array([.1,.1,30.])[:,None],velocity)
    for frame_name,frame,quadrature,steps in [('moving',boost,3,n) for n in (25,50,100)]+[('lab',np.zeros(3),3,100),('moving',boost,2,100)]:
        source=FixedFrameAtomicMoments(network,reference.neutral_mass,frame);model=MaterialReactionOverlap(source,density,quadrature=quadrature)
        state=model.from_material(material,initial_material,[0,1,2]);initial=state;times=[0.];states=[state.moments.copy()]
        budgets=[model.totals(state)];initial_fields=model.fields(state);failure=None;moment_error=0.;population_error=0.;energy_error=0.
        charge_error=0.;current_error=0.;maximum_drift=0.
        for step in range(1,steps+1):
            try:candidate=model.step(state,duration/steps)
            except (ValueError,RuntimeError) as error:
                failure=dict(step=step,time_s=times[-1],message=str(error));break
            state=candidate;times.append(step*duration/steps);states.append(state.moments.copy());budgets.append(model.totals(state))
            exact=exact_moments(source,initial.moments,times[-1],density);delta=state.moments-exact
            population_error=max(population_error,float(abs(delta[...,model.number]).max()))
            moment_error=max(moment_error,float(abs(delta[...,source.momentum]).max()))
            energy_error=max(energy_error,float(abs(delta[...,source.energy]).max()))
            fields=model.fields(state);charge_error=max(charge_error,float(abs(fields['charge_c_m3']-initial_fields['charge_c_m3']).max()))
            current_error=max(current_error,float(np.linalg.norm(fields['current_a_m2']-initial_fields['current_a_m2'],axis=1).max()))
        name=f'{frame_name}_q{quadrature}_t{steps}';budget={key:np.asarray([row[key] for row in budgets]) for key in budgets[0]}
        final=model.fields(state);x=state.points.ravel();born=(abs(x)>.005)&(abs(x)<.010);empty=abs(x)>.015;unseeded=(abs(x)>.010)&(abs(x)<.015)
        record=dict(name=name,success=failure is None,failure=failure,steps=len(times)-1,requested_steps=steps,completed_time_s=times[-1],
            maximum_analytic_population_error=population_error,maximum_analytic_scaled_momentum_error=moment_error,
            maximum_analytic_translational_energy_error_ev=energy_error,maximum_local_charge_change_c_m3=charge_error,
            maximum_current_change_a_m2=current_error,
            maximum_mass_change_kg_m2=float(abs(budget['mass_kg_m2']-budget['mass_kg_m2'][0]).max()),
            maximum_nucleus_change_m2=float(abs(budget['nuclei_m2']-budget['nuclei_m2'][0]).max()),
            maximum_momentum_change_kg_per_m_s=float(np.linalg.norm(budget['momentum_kg_per_m_s']-budget['momentum_kg_per_m_s'][0],axis=1).max()),
            maximum_matter_energy_change_j_m2=float(abs(budget['matter_binding_photon_energy_j_m2']-budget['matter_binding_photon_energy_j_m2'][0]).max()),
            maximum_electric_change_v_m=float(abs(budget['electric_edges_v_m']-budget['electric_edges_v_m'][0]).max()),
            born_ions_outside_old_support_m2=float(np.sum(state.weights.ravel()[born]*final['number_m3'][born,1])),
            maximum_exterior_number_m3=float(final['number_m3'][empty].max(initial=0.)),
            maximum_unseeded_electron_number_m3=float(final['number_m3'][unseeded,2].max(initial=0.)),
            final_electrons_m2=float(budget['numbers_m2'][-1,-1]),frame_velocity_m_s=frame.tolist(),quadrature=quadrature)
        np.savez_compressed(out/(name+'.npz'),time_s=times,moments=states,points_m=state.points,weights_m=state.weights,edges_m=state.edges,
            initial_material_faces_m=faces,initial_material_numbers_m2=number,initial_material_velocity_m_s=velocity,
            final_number_m3=final['number_m3'],final_current_a_m2=final['current_a_m2'],**budget)
        report['cases'].append(record);histories[name]=(source,model,initial,state,budget);save();print(record,flush=True)
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    source,model,initial,state,budget=histories['moving_q3_t100'];x=state.points.ravel();first=model.fields(initial);last=model.fields(state)
    for k,label in enumerate(('Neutral','Ion','Electron')):
        axes[0].plot(x*1000,first['number_m3'][:,k]/density,'--',label=f'{label}: initial')
        axes[0].plot(x*1000,last['number_m3'][:,k]/density,label=f'{label}: reacted')
    axes[1].plot(x*1000,first['current_a_m2'][:,1],label='Initial');axes[1].plot(x*1000,last['current_a_m2'][:,1],'--',label='Reacted')
    axes[2].plot(np.linspace(0,duration,len(budget['energy_j_m2']))*1e6,budget['matter_binding_photon_energy_j_m2']-budget['matter_binding_photon_energy_j_m2'][0])
    axes[0].set(xlabel='Position [mm]',ylabel='Number density / reference',title='Actual product support grows')
    axes[1].set(xlabel='Position [mm]',ylabel='Transverse current [A/m²]',title='Equal/opposite ionization birth current')
    axes[2].set(xlabel='Reaction-source time [µs]',ylabel='Matter + binding energy change [J/m²]',title='Independent laboratory energy budget')
    for ax in axes:ax.grid(alpha=.2)
    for ax in axes[:2]:ax.legend(fontsize=7)
    fig.savefig(out/'reaction_overlap.png',dpi=170);fig.savefig(out/'reaction_overlap.pdf');plt.close(fig)
    if identity!=provenance(['examples.material_reactions']):raise RuntimeError('Source changed during reaction-overlap run')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');args=parser.parse_args()
    result=run(args.destination);raise SystemExit(0 if result['success'] else 1)
