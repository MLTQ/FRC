"""Actual moving reaction limits and unequal-mesh electron-ion dynamics."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass,electron_mass,elementary_charge as EV
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.species_material import MaterialSpecies
from plasma_bench.material_species_set import MaterialSpeciesSet,MaterialSetState
from plasma_bench.material_reactions import MaterialReactionOverlap,ReactionOverlapState
from plasma_bench.material_reaction_projection import project_reactions
from plasma_bench.reactive_material import ReactiveMaterial,MaterialReservoirs
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance,freeze


class NeutralDecayLimit:
    """Declared two-level spontaneous-decay limit; no fitted atomic-data claim."""
    data={'charge':0,'atomic_number':1}
    energy_ev=np.array([0.,3.])
    radiation=SimpleNamespace(lower=np.array([0]),upper=np.array([1]),A=np.array([2e7]),photon_energy=np.array([3.])*EV)

    def generator(self,temperature,electron_density,*,radiation=False):
        if electron_density:raise ValueError('This neutral mathematical limit has no supplied electron excitation rates')
        return np.zeros((4,4))


def recombination():
    return RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',A_cm3_s=1e-12,
        reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))


def write_states(path,models,states):
    arrays={}
    for k,(model,state) in enumerate(zip(models,states)):
        for b,(m,s) in enumerate(zip(model.models,state.blocks)):
            for name,value in [('faces_m',s.faces),('relative_velocity_m_s',s.velocity),('thermal_j_m2',s.thermal),('numbers_m2',m.numbers)]:
                arrays[f'step{k}_block{b}_{name}']=value
            if hasattr(m,'velocity_indices'):
                arrays.pop(f'step{k}_block{b}_numbers_m2')
                arrays[f'step{k}_block{b}_numbers_divided_by_number_scale']=m.relative_numbers
                arrays[f'step{k}_block{b}_thermal_divided_by_number_scale']=arrays.pop(f'step{k}_block{b}_thermal_j_m2')
                arrays[f'step{k}_block{b}_number_scale_m2']=np.asarray(m.number_scale)
                arrays[f'step{k}_block{b}_velocity_indices']=m.velocity_indices
                arrays[f'step{k}_block{b}_reference_faces_m']=m.reference_faces
        arrays[f'step{k}_species_indices']=model.species_indices
        arrays[f'step{k}_frame_velocity_m_s']=model.reference_velocity
    np.savez_compressed(path,**arrays)


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh moving reactive output required')
    identity=provenance(['examples.reactive_material']);out.mkdir(parents=True);freeze(out,identity);start=time.monotonic()
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],
        scope='Complete moving cold neutral radiative limit with variable level masses and photons; unequal-mesh physical-mass20ns electron-ion transport limit; not general hot charged reacting formation adoption.')
    def save():
        report['elapsed_seconds']=time.monotonic()-start;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    network=AtomicNetwork(ConstantIonization(),recombination(),{0:NeutralDecayLimit()},ionization_enabled=False,recombination_enabled=False)
    frame=np.array([1000.,-2000.,500.]);source=FixedFrameAtomicMoments(network,1.00794*atomic_mass,frame)
    for steps in (20,40,80):
        faces=np.linspace(-.02,.02,9)[None,:];single=MaterialSpecies(1e19*np.diff(faces),[source.mass[1]],[0.],(-.1,.1))
        model=MaterialSpeciesSet([single],[1],reference_velocity=frame);state=MaterialSetState((single.state(faces,0.),))
        solver=ReactiveMaterial(source,1e19);reservoirs=MaterialReservoirs();duration=40e-9;times=[0.];models=[model];states=[state];ledger=[solver.totals(model,state,reservoirs)];failure=None
        for k in range(1,steps+1):
            try:model,state,reservoirs=solver.step(model,state,reservoirs,duration/steps)
            except (ValueError,RuntimeError) as error:
                failure=dict(step=k,time_s=times[-1],message=str(error))
                if solver.last_overlap is not None:
                    r=solver.last_overlap;np.savez_compressed(out/f'decay_t{steps}_rejected_overlap.npz',edges_m=r.edges,points_m=r.points,weights_m=r.weights,moments=r.moments)
                break
            times.append(k*duration/steps);models.append(model);states.append(state);ledger.append(solver.totals(model,state,reservoirs))
        number=np.asarray([row['atomic_numbers_m2'] for row in ledger]);energy=np.array([r['full_energy_j_m2'] for r in ledger]);photons=np.array([r['photon_energy_j_m2'] for r in ledger])
        initial_number=4e17;exact_upper=initial_number*np.exp(-2e7*np.asarray(times));position_error=0.;velocity_error=0.
        for t,m,s in zip(times,models,states):
            for block in s.blocks:
                position_error=max(position_error,float(abs(block.faces-(faces+frame[2]*t)).max()))
                velocity_error=max(velocity_error,float(abs(block.velocity).max()))
        record=dict(name=f'decay_t{steps}',success=failure is None,failure=failure,completed_time_s=times[-1],requested_time_s=duration,
            maximum_relative_population_error=float(abs(number[:,1]-exact_upper).max()/initial_number),
            maximum_position_error_m=position_error,maximum_velocity_error_m_s=velocity_error,
            maximum_energy_change_j_m2=float(abs(energy-energy[0]).max()),maximum_photon_error_j_m2=float(abs(photons-3*EV*(initial_number-exact_upper)).max()),
            maximum_electron_number_m2=float(number[:,-1].max()),maximum_cells=sum(b.cells for b in model.models))
        write_states(out/(record['name']+'.npz'),models,states)
        np.savez_compressed(out/(record['name']+'_budgets.npz'),time_s=times,atomic_numbers_m2=number,full_energy_j_m2=energy,photon_energy_j_m2=photons)
        report['cases'].append(record);save();print(record,flush=True)
        axes[0].plot(np.asarray(times)*1e9,number[:,1]/initial_number,label=f'{steps} steps')
    t=np.linspace(0,40e-9,200);axes[0].plot(t*1e9,np.exp(-2e7*t),'k--',label='Analytic decay + translation')
    finals=[]
    for cells in (32,64,128):
        models=[];blocks=[];masses=[39.948*atomic_mass-electron_mass,electron_mass]
        for n,mass,charge,temp in [(cells,masses[0],1.,.2),(2*cells,masses[1],-1.,2.)]:
            faces=-.01*np.cos(np.linspace(0,np.pi,n+1))[None,:]
            m=MaterialSpecies(1e14*np.diff(faces),[mass],[charge],(-.2,.2),maximum_plasma_phase=.04*32/(2*cells),numerical_viscosity=(.5,2.),viscosity_limiter='minmod')
            models.append(m);blocks.append(m.state(faces,temp))
        model=MaterialSpeciesSet(models,[0,1]);state=MaterialSetState(tuple(blocks));times=[0.];states=[state];energies=[model.totals(state)['energy_j_m2']];failure=None;duration=20e-9
        while times[-1]<duration:
            dt=min(model.timestep(state),duration-times[-1])
            try:candidate=model.step(state,dt)
            except (ValueError,RuntimeError) as error:failure=dict(time_s=times[-1],dt_s=dt,message=str(error));break
            state=candidate;times.append(min(times[-1]+dt,duration));states.append(state);energies.append(model.totals(state)['energy_j_m2'])
        record=dict(name=f'plasma_i{cells}_e{2*cells}',success=failure is None,failure=failure,completed_time_s=times[-1],requested_time_s=duration,steps=len(times)-1,
            maximum_relative_energy_error=float(np.max(abs(np.asarray(energies)-energies[0]))/energies[0]),fronts_m=[s.faces[0,[0,-1]].tolist() for s in state.blocks])
        write_states(out/(record['name']+'.npz'),[model]*len(states),states);np.savez_compressed(out/(record['name']+'_budgets.npz'),time_s=times,energy_j_m2=energies)
        report['cases'].append(record);save();print(record,flush=True);finals.append((model,state))
        x=np.linspace(-.07,.07,8001);f=model.sample(state,x);axes[1].plot(x*1000,f['electric'],label=f'{cells}/{2*cells} cells');axes[2].plot(x*1000,f['current'][:,2],label=f'{cells}/{2*cells} cells')
    comparisons=[]
    for (m,s),(n,t) in zip(finals[:-1],finals[1:]):
        edges=np.unique(np.r_[m.bounds,s.faces,t.faces]);nodes=np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
        x=(edges[:-1,None]+np.diff(edges)[:,None]*nodes).ravel();w=np.repeat(np.diff(edges)/2,2);a=m.sample(s,x);b=n.sample(t,x)
        comparisons.append(dict(coarse_counts=[q.cells for q in m.models],fine_counts=[q.cells for q in n.models],
            electric_relative_l2=float(np.sqrt(np.sum(w*(a['electric']-b['electric'])**2)/np.sum(w*b['electric']**2))),
            current_relative_l2=float(np.sqrt(np.sum(w[:,None]*(a['current']-b['current'])**2)/np.sum(w[:,None]*b['current']**2)))))
    report['plasma_spatial_comparisons']=comparisons
    for ax,title in zip(axes,('Moving reactive upper population / initial','Electric field [V/m]','Current [A/m²]')):ax.set_ylabel(title);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[0].set_xlabel('Time [ns]')
    for ax in axes[1:]:ax.set_xlabel('Position [mm]')
    fig.savefig(out/'moving_reactions.png',dpi=170);fig.savefig(out/'moving_reactions.pdf');plt.close(fig)
    if identity!=provenance(['examples.reactive_material']):raise RuntimeError('Source changed during actual reactive-material run')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');args=parser.parse_args()
    result=run(args.destination);raise SystemExit(0 if result['success'] else 1)
