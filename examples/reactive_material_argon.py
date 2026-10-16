"""Attempt complete hot charged reacting argon motion without hidden repairs."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.species_material import MaterialSpecies
from plasma_bench.material_subcells import MaterialSubcells
from plasma_bench.material_species_set import MaterialSpeciesSet,MaterialSetState
from plasma_bench.reactive_material import ReactiveMaterial,MaterialReservoirs
from .reactive_material import write_states
from .implicit_plasma_records import provenance,freeze


def run(destination,*,separate_velocity=False,density_width_fraction=0.):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh charged reacting output required')
    if not np.isfinite(density_width_fraction) or density_width_fraction<0:raise ValueError('Nonnegative density width fraction required')
    data=Path(__file__).resolve().parents[1]/'data/atomic'
    ionization=ImpactIonization.from_file(data/'argon/argon.json');rr=RadiativeRecombination.from_file(data/'recombination/argon_rr.json')
    network=AtomicNetwork(ionization,rr);source=FixedFrameAtomicMoments(network,39.948*atomic_mass,[0.,0.,0.])
    identity=provenance(['examples.reactive_material_argon']);out.mkdir(parents=True);freeze(out,identity);start=time.monotonic()
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],maximum_occupied_cells=10000,maximum_steps=5000,minimum_timestep_s=2e-17,
        separate_velocity=separate_velocity,density_width_fraction=density_width_fraction,geometry_precision_bits=np.finfo(np.longdouble).nmant+1 if separate_velocity else 53,
        inputs=dict(ionization_sha256=ionization.input_sha256,recombination_sha256=rr.input_sha256),
        scope='Physical argon masses, existing ground charge ladder/DI/RR, neutral+ion+electron finite supports and common electrostatic field. Full20ns requested at32/64/128. Real rate/model-range limits remain; no collisions/excitation/kinetic validation or formation adoption.')
    def save():
        report['elapsed_seconds']=time.monotonic()-start;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    for cells in (32,64,128):
        faces=-.01*np.cos(np.linspace(0,np.pi,cells+1))[None,:];ids=[network.ground[0],network.ground[1],network.count];models=[];blocks=[]
        for index,temp in zip(ids,[.1,.2,2.]):
            controls=dict(maximum_plasma_phase=.04*32/cells,numerical_viscosity=(.5,2.),viscosity_limiter='minmod')
            if separate_velocity:
                m=MaterialSubcells(1e14*np.diff(faces)[0],source.mass[index],np.r_[network.charge,-1.][index],(-.2,.2),faces[0],cells,**controls)
            else:m=MaterialSpecies(1e14*np.diff(faces),[source.mass[index]],[np.r_[network.charge,-1.][index]],(-.2,.2),**controls)
            models.append(m);blocks.append(m.state(faces,temp))
        model=MaterialSpeciesSet(models,ids);state=MaterialSetState(tuple(blocks));ledger=MaterialReservoirs()
        minimum_width=density_width_fraction*.01*(1-np.cos(np.pi/cells))
        solver=ReactiveMaterial(source,1e14,maximum_plasma_phase=.04*32/cells,coalesce_ulps=8,minimum_density_width=minimum_width,**(dict(velocity_cells=cells) if separate_velocity else {}))
        times=[0.];models_history=[model];states=[state];budgets=[solver.totals(model,state,ledger)];failure=None;duration=20e-9
        while times[-1]<duration:
            occupied=sum(m.cells for m in model.models)
            if occupied>report['maximum_occupied_cells'] or len(times)-1>=report['maximum_steps']:
                failure=dict(time_s=times[-1],message='Declared computational limit reached; requested physical duration not completed',occupied_cells=occupied);break
            dt=min(.5*model.timestep(state),duration-times[-1])
            if not np.isfinite(dt) or dt<=0 or times[-1]+dt==times[-1]:
                failure=dict(time_s=times[-1],dt_s=float(dt),message='Actual material timestep is nonfinite or cannot advance physical time');break
            if dt<report['minimum_timestep_s']:
                failure=dict(time_s=times[-1],dt_s=float(dt),message='Declared computational timestep limit reached; requested20ns not completed',occupied_cells=occupied);break
            try:new_model,new_state,new_ledger=solver.step(model,state,ledger,dt)
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
                failure=dict(time_s=times[-1],dt_s=dt,message=str(error),occupied_cells=occupied)
                if solver.last_overlap is not None:
                    r=solver.last_overlap;np.savez_compressed(out/f'argon_n{cells}_rejected_overlap.npz',edges_m=r.edges,points_m=r.points,weights_m=r.weights,moments=r.moments)
                break
            model,state,ledger=new_model,new_state,new_ledger;times.append(min(times[-1]+dt,duration));models_history.append(model);states.append(state);budgets.append(solver.totals(model,state,ledger))
            if len(times)==2 or (len(times)-1)%20==0:print('argon',cells,'step',len(times)-1,'time_s',times[-1],'cells',sum(m.cells for m in model.models),flush=True)
        energy=np.array([r['full_energy_j_m2'] for r in budgets]);number=np.stack([r['atomic_numbers_m2'] for r in budgets])
        record=dict(name=f'argon_n{cells}',success=failure is None,failure=failure,steps=len(times)-1,completed_time_s=times[-1],requested_time_s=duration,
            minimum_density_width_m=minimum_width,observed_projection_field_change_j_m2=ledger.projection_field_change_j_m2,
            final_occupied_cells=sum(m.cells for m in model.models),maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
            final_electrons_m2=float(number[-1,-1]),initial_electrons_m2=float(number[0,-1]))
        widths=[np.diff(s.faces[0]) for s in state.blocks];b=int(np.argmin([w.min() for w in widths]));j=int(widths[b].argmin())
        record['smallest_cell']=dict(species=int(model.species_indices[b]),block=b,cell=j,width_m=float(widths[b][j]),
            endpoints_m=state.blocks[b].faces[0,j:j+2].astype(float).tolist(),number_m2=float(model.models[b].numbers[0,j]),
            support_m=state.blocks[b].faces[0,[0,-1]].astype(float).tolist(),minimum_thermal_j_m2=float(state.blocks[b].thermal.min()*getattr(model.models[b],'number_scale',1.)))
        record['occupied_species']=model.species_indices.tolist()
        write_states(out/(record['name']+'.npz'),models_history,states)
        np.savez_compressed(out/(record['name']+'_budgets.npz'),time_s=times,atomic_numbers_m2=number,full_energy_j_m2=energy,
                            observed_projection_field_change_j_m2=[r['observed_projection_field_change_j_m2'] for r in budgets],
                            photon_energy_j_m2=[r['photon_energy_j_m2'] for r in budgets])
        report['cases'].append(record);save();print(record,flush=True)
    if identity!=provenance(['examples.reactive_material_argon']):raise RuntimeError('Source changed during charged reactive attempt')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--separate-velocity',action='store_true');parser.add_argument('--density-width-fraction',type=float,default=0.);args=parser.parse_args()
    result=run(args.destination,separate_velocity=args.separate_velocity,density_width_fraction=args.density_width_fraction);raise SystemExit(0 if result['success'] else 1)
