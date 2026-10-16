"""Complete independent physical limits for the two-space material solver."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass,electron_mass as ME,elementary_charge as EV,epsilon_0 as EPS0
from plasma_bench.material_subcells import MaterialSubcells
from plasma_bench.material_species_set import MaterialSpeciesSet,MaterialSetState
from plasma_bench.reactive_material import ReactiveMaterial,MaterialReservoirs
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from .reactive_material import NeutralDecayLimit,recombination,write_states
from .atomic_moments import ConstantIonization
from .atomic_moving import charged_expansion_time
from .implicit_plasma_records import provenance,freeze


def subdivided_faces(cells,a):
    coarse=-a*np.cos(np.linspace(0,np.pi,cells+1))
    return np.r_[(coarse[:-1,None]+np.diff(coarse)[:,None]*np.arange(4)/4).ravel(),coarse[-1]]


def save_fixed(path,model,states,times,budgets):
    arrays=dict(time_s=times,energy_j_m2=[r['energy_j_m2'] for r in budgets],
        momentum_kg_per_m_s=[r['momentum_kg_per_m_s'] for r in budgets],species_indices=model.species_indices)
    for b,m in enumerate(model.models):
        prefix=f'block{b}_'
        for key,value in dict(numbers_divided_by_number_scale=m.relative_numbers,number_scale_m2=m.number_scale,
            reference_faces_m=m.reference_faces,velocity_indices=m.velocity_indices,
            faces_m=[s.blocks[b].faces for s in states],velocity_m_s=[s.blocks[b].velocity for s in states],
            thermal_divided_by_number_scale=[s.blocks[b].thermal for s in states]).items():arrays[prefix+key]=np.asarray(value)
    np.savez_compressed(path,**arrays)


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh subcell limit output required')
    identity=provenance(['examples.material_subcells']);out.mkdir(parents=True);freeze(out,identity);began=time.monotonic()
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],geometry_precision_bits=np.finfo(np.longdouble).nmant+1,
        scope='Full moving neutral radiative decay, independently soluble warm charged affine expansion, and full20ns ion/electron expansion. Distinct fine density/coarse velocity spaces; no claim of successful hot reacting argon or formation.')
    def save():
        report['elapsed_seconds']=time.monotonic()-began;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    network=AtomicNetwork(ConstantIonization(),recombination(),{0:NeutralDecayLimit()},ionization_enabled=False,recombination_enabled=False)
    frame=np.array([1000.,-2000.,500.]);source=FixedFrameAtomicMoments(network,1.00794*atomic_mass,frame)
    for steps in (20,40,80):
        faces=subdivided_faces(8,.02);block=MaterialSubcells(1e19*np.diff(faces),source.mass[1],0.,(-.1,.1),faces,8)
        model=MaterialSpeciesSet([block],[1],reference_velocity=frame);state=MaterialSetState((block.state(faces[None,:],0.),))
        solver=ReactiveMaterial(source,1e19,velocity_cells=8,coalesce_ulps=8);ledger=MaterialReservoirs();duration=40e-9
        times=[0.];models=[model];states=[state];budgets=[solver.totals(model,state,ledger)];failure=None
        for k in range(1,steps+1):
            try:model,state,ledger=solver.step(model,state,ledger,duration/steps)
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:failure=str(error);break
            times.append(k*duration/steps);models.append(model);states.append(state);budgets.append(solver.totals(model,state,ledger))
        number=np.array([r['atomic_numbers_m2'] for r in budgets]);energy=np.array([r['full_energy_j_m2'] for r in budgets]);photons=np.array([r['photon_energy_j_m2'] for r in budgets]);exact=4e17*np.exp(-2e7*np.asarray(times))
        position=max(float(abs(s.blocks[0].faces[0,[0,-1]]-np.array([-.02,.02])-frame[2]*t).max()) for s,t in zip(states,times))
        row=dict(name=f'decay_t{steps}',success=failure is None,failure=failure,requested_time_s=duration,completed_time_s=times[-1],steps=len(times)-1,
            maximum_relative_population_error=float(abs(number[:,1]-exact).max()/4e17),maximum_front_error_m=position,
            maximum_energy_change_j_m2=float(abs(energy-energy[0]).max()),maximum_photon_error_j_m2=float(abs(photons-3*EV*(4e17-exact)).max()))
        write_states(out/(row['name']+'.npz'),models,states);np.savez_compressed(out/(row['name']+'_budgets.npz'),time_s=times,atomic_numbers_m2=number,energy_j_m2=energy,photon_energy_j_m2=photons)
        report['cases'].append(row);save();print(row,flush=True);axes[0].plot(np.asarray(times)*1e9,number[:,1]/4e17,label=f'{steps} steps')
    x=np.linspace(0,40e-9,200);axes[0].plot(x*1e9,np.exp(-2e7*x),'k--',label='Analytic decay')
    finals=[];mass=39.948*atomic_mass-ME
    for kind in ('charged','plasma'):
        for cells in (32,64,128):
            a=.025 if kind=='charged' else .01;density=1e12 if kind=='charged' else 1e14;faces=subdivided_faces(cells,a)
            params=[(mass,1.,1.)] if kind=='charged' else [(mass,1.,.2),(ME,-1.,2.)];blocks=[];models=[]
            for m,q,temp in params:
                model=MaterialSubcells(density*np.diff(faces),m,q,(-.2,.2),faces,cells,
                    left_electric=-EV*density*a/EPS0 if kind=='charged' else 0.,maximum_plasma_phase=.04*32/cells,numerical_viscosity=(.5,2.),viscosity_limiter='minmod')
                T=temp*(1-(faces[:-1]**2+faces[:-1]*faces[1:]+faces[1:]**2)/(3*a*a)) if kind=='charged' else temp
                models.append(model);blocks.append(model.state(faces[None,:],T))
            model=MaterialSpeciesSet(models,range(len(models)));state=MaterialSetState(tuple(blocks));states=[state];times=[0.];budgets=[model.totals(state)];failure=None
            duration,scale_rate=charged_expansion_time(1.5,2*EV/(mass*a*a),EV**2*density/(EPS0*mass)) if kind=='charged' else (20e-9,None)
            while times[-1]<duration:
                dt=min(model.timestep(state),duration-times[-1])
                try:state=model.step(state,dt)
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:failure=str(error);break
                times.append(min(times[-1]+dt,duration));states.append(state);budgets.append(model.totals(state))
            energy=np.array([r['energy_j_m2'] for r in budgets]);momentum=np.array([r['momentum_kg_per_m_s'] for r in budgets])
            row=dict(name=f'{kind}_n{cells}',success=failure is None,failure=failure,velocity_cells=cells,density_cells=4*cells,
                requested_time_s=duration,completed_time_s=times[-1],steps=len(times)-1,
                maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),maximum_momentum_change=float(abs(momentum-momentum[0]).max()))
            if kind=='charged' and failure is None:
                row.update(maximum_position_error_m=float(abs(state.blocks[0].faces[0]-1.5*faces).max()),
                    maximum_velocity_error_m_s=float(abs(state.blocks[0].velocity[0,:,2]-scale_rate*faces[model.models[0].velocity_indices]).max()))
            save_fixed(out/(row['name']+'.npz'),model,states,times,budgets);report['cases'].append(row);save();print(row,flush=True)
            if kind=='plasma':
                finals.append((model,state,times[-1]));x=np.linspace(-.07,.07,4001);f=model.sample(state,x)
                label=f'{cells} velocity cells'+(f' STOPPED at{times[-1]*1e9:.4g}ns' if failure else '')
                axes[1].plot(x*1000,f['electric'],label=label);axes[2].plot(x*1000,f['current'][:,2],label=label)
    report['plasma_spatial_comparisons']=[]
    for (m,s,ta),(n,t,tb) in zip(finals[:-1],finals[1:]):
        if ta!=20e-9 or tb!=20e-9:continue
        edges=np.unique(np.r_[m.bounds,s.faces,t.faces]);nodes=np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
        x=(edges[:-1,None]+np.diff(edges)[:,None]*nodes).ravel();w=np.repeat(np.diff(edges)/2,2);a=m.sample(s,x);b=n.sample(t,x)
        report['plasma_spatial_comparisons'].append(dict(coarse_velocity_cells=m.models[0].velocity_cells,fine_velocity_cells=n.models[0].velocity_cells,
            electric_relative_l2=float(np.sqrt(np.sum(w*(a['electric']-b['electric'])**2)/np.sum(w*b['electric']**2))),
            current_relative_l2=float(np.sqrt(np.sum(w[:,None]*(a['current']-b['current'])**2)/np.sum(w[:,None]*b['current']**2)))))
    for ax,label in zip(axes,('Reactive upper population / initial','Electric field [V/m]','Current [A/m²]')):ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[0].set_xlabel('Time [ns]')
    for ax in axes[1:]:ax.set_xlabel('Position [mm]')
    fig.savefig(out/'subcell_limits.png',dpi=170);fig.savefig(out/'subcell_limits.pdf');plt.close(fig)
    if provenance(['examples.material_subcells'])!=identity:raise RuntimeError('Source changed during complete subcell limits')
    report['terminal']=True;report['success']=all(r['success'] for r in report['cases']);save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');a=p.parse_args();r=run(a.destination);raise SystemExit(0 if r['success'] else 1)
