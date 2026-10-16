"""Full warm affine motion with neutral level decay and density coarsening."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass,elementary_charge as EV
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.material_subcells import MaterialSubcells
from plasma_bench.material_species_set import MaterialSpeciesSet,MaterialSetState
from plasma_bench.reactive_material import ReactiveMaterial,MaterialReservoirs
from .reactive_material import NeutralDecayLimit,recombination,write_states
from .material_subcells import subdivided_faces
from .atomic_moments import ConstantIonization
from .atomic_moving import charged_expansion_time
from .implicit_plasma_records import provenance,freeze


class WarmDecay(NeutralDecayLimit):
    radiation=SimpleNamespace(lower=np.array([0]),upper=np.array([1]),A=np.array([2e5]),photon_energy=np.array([3.])*EV)


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh warm reacting affine output required')
    identity=provenance(['examples.reactive_material_affine']);out.mkdir(parents=True);freeze(out,identity);began=time.monotonic()
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],
        scope='Independent full warm neutral affine expansion with actual radiative level transfer; declared mathematical two-level hydrogen-mass atom, not argon experimental or formation validation.')
    def save():
        report['elapsed_seconds']=time.monotonic()-began;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();network=AtomicNetwork(ConstantIonization(),recombination(),{0:WarmDecay()},ionization_enabled=False,recombination_enabled=False)
    mass=1.00794*atomic_mass;source=FixedFrameAtomicMoments(network,mass,[0.,0.,0.]);a=.025;density=1e19;temperature=.1;initial_number=2*a*density
    duration,rate=charged_expansion_time(1.5,2*EV*temperature/(mass*a*a),0.)
    report.update(duration_s=duration,analytic_final_scale=1.5,analytic_final_scale_rate_s=rate,decay_rate_s=2e5)
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for cells,fraction in [(16,.1),(32,.1),(64,.1),(32,.05),(32,.025)]:
        faces=subdivided_faces(cells,a);T=temperature*(1-(faces[:-1]**2+faces[:-1]*faces[1:]+faces[1:]**2)/(3*a*a))
        m=MaterialSubcells(density*np.diff(faces),mass,0.,(-.2,.2),faces,cells,numerical_viscosity=(.5,2.),viscosity_limiter='minmod')
        model=MaterialSpeciesSet([m],[1]);state=MaterialSetState((m.state(faces[None,:],T),));ledger=MaterialReservoirs()
        width=fraction*a*(1-np.cos(np.pi/cells));solver=ReactiveMaterial(source,density,velocity_cells=cells,coalesce_ulps=8,minimum_density_width=width)
        times=[0.];models=[model];states=[state];budgets=[solver.totals(model,state,ledger)];failure=None;maximum_cells=m.cells
        while times[-1]<duration:
            dt=min(.5*model.timestep(state),duration-times[-1])
            if sum(q.cells for q in model.models)>10000 or len(times)>10000:
                failure=dict(time_s=times[-1],message='Declared10000-cell/step computational bound');break
            try:new_model,new_state,new_ledger=solver.step(model,state,ledger,dt)
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
                failure=dict(time_s=times[-1],dt_s=dt,message=str(error));break
            model,state,ledger=new_model,new_state,new_ledger;times.append(min(times[-1]+dt,duration));models.append(model);states.append(state);budgets.append(solver.totals(model,state,ledger))
            maximum_cells=max(maximum_cells,sum(q.cells for q in model.models))
        number=np.asarray([r['atomic_numbers_m2'] for r in budgets]);energy=np.asarray([r['full_energy_j_m2'] for r in budgets]);photons=np.asarray([r['photon_energy_j_m2'] for r in budgets]);exact=initial_number*np.exp(-2e5*np.asarray(times))
        row=dict(name=f'affine_n{cells}_f{fraction:g}',velocity_cells=cells,density_width_fraction=fraction,minimum_density_width_m=width,
            success=failure is None,failure=failure,steps=len(times)-1,completed_time_s=times[-1],requested_time_s=duration,maximum_occupied_cells=maximum_cells,
            maximum_relative_population_error=float(abs(number[:,1]-exact).max()/initial_number),maximum_energy_change_j_m2=float(abs(energy-energy[0]).max()),
            maximum_photon_error_j_m2=float(abs(photons-3*EV*(initial_number-exact)).max()),observed_projection_field_change_j_m2=ledger.projection_field_change_j_m2)
        if failure is None:
            row['maximum_front_error_m']=max(float(abs(s.faces[0,[0,-1]]-1.5*np.array([-a,a])).max()) for s in state.blocks)
            row['maximum_velocity_error_m_s']=max(float(abs(s.velocity[0,:,2]-rate*s.faces[0,m.velocity_indices]/1.5).max()) for m,s in zip(model.models,state.blocks))
        else:
            if solver.last_overlap is not None:
                r=solver.last_overlap;np.savez_compressed(out/(row['name']+'_failed_overlap.npz'),edges_m=r.edges,points_m=r.points,weights_m=r.weights,moments=r.moments)
        write_states(out/(row['name']+'.npz'),models,states);np.savez_compressed(out/(row['name']+'_budgets.npz'),time_s=times,atomic_numbers_m2=number,full_energy_j_m2=energy,photon_energy_j_m2=photons)
        report['cases'].append(row);save();print(row,flush=True)
        label=f'{cells} velocity cells, f={fraction:g}'+(' STOPPED' if failure else '')
        axes[0].plot(np.asarray(times)*1e6,number[:,1]/initial_number,label=label)
        axes[1].plot(np.asarray(times)*1e6,(energy-energy[0])/energy[0],label=label)
    t=np.linspace(0,duration,200);axes[0].plot(t*1e6,np.exp(-2e5*t),'k--',label='Analytic decay')
    for ax in axes:ax.set_xlabel('Time [µs]');ax.grid(alpha=.2);ax.legend(fontsize=6)
    axes[0].set_ylabel('Excited level / initial population');axes[1].set_ylabel('Physical total-energy change / initial')
    fig.savefig(out/'warm_reacting_affine.png',dpi=170);fig.savefig(out/'warm_reacting_affine.pdf');plt.close(fig)
    if identity!=provenance(['examples.reactive_material_affine']):raise RuntimeError('Source changed during warm reacting affine run')
    report['terminal']=True;report['success']=all(r['success'] for r in report['cases']);save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');a=p.parse_args();r=run(a.destination);raise SystemExit(0 if r['success'] else 1)
