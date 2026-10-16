"""Continue an actual constant-current pulse interval to localize energy loss."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from plasma_bench.material_experiment import build
from plasma_bench.material_plasma_samples import material_samples
from .material_coil_plot import state_at
from .implicit_plasma_records import provenance, freeze, sha


def arrays(state):
    return dict(magnetic=state.magnetic, **{k:getattr(state.fluid,k) for k in state.fluid.__dataclass_fields__})


def run(source, destination, first=1450, last=1500, *, central_integrator=None, factors=(1,2,4)):
    source=Path(source);out=Path(destination)
    if out.exists():raise ValueError('Fresh output directory required')
    config=json.loads((source/'config.json').read_text())
    if central_integrator is not None:config['central_integrator']=central_integrator
    if not factors or sorted(set(factors))!=list(factors) or min(factors)<1:
        raise ValueError('Increasing positive integer refinement factors required')
    original_dt=config['evolution']['duration_s']/config['evolution']['steps']
    start=first*original_dt;end=last*original_dt
    if first>=last or config['integrator']!='energy_split':raise ValueError('Ordered saved energy-split interval required')
    initial=state_at(source,first);saved_final=state_at(source,last)
    identity=provenance(['examples.material_hold_interval']);out.mkdir(parents=True);freeze(out,identity)
    shutil.copy2(source/f'step_{first:05}.npz',out/'initial.npz')
    shutil.copy2(source/f'step_{last:05}.npz',out/'original_final.npz')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    began=time.monotonic();records=[];endpoints=[]
    report=dict(terminal=False,source_sha256=identity['source_sha256'],source_run=str(source.resolve()),
        source_initial_sha256=sha(out/'initial.npz'),source_final_sha256=sha(out/'original_final.npz'),
        start_s=start,end_s=end,scope='Actual saved-state constant-current hold continuations; independent invariant is constant total energy, with unchanged physical equations and exact vacuum. This local interval does not replace full-pulse or spatial convergence.',
        central_integrator=config.get('central_integrator','rk4'),factors=list(factors),directional_power=[],runs=records)
    def save():
        report['elapsed_seconds']=time.monotonic()-began
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    _,probe,_=build(config)
    # A directional difference of actual total energy, independent of the coded power pairing.
    for step in (first,(first+last)//2,last):
        state=state_at(source,step);t=step*original_dt;a=probe.response(state,time=t)
        derivative=probe.derivative(state,time=t)
        row=dict(step=step,time_s=t,external_power_w=a['external_power'],pair_power_residual_w=a['pair_power_residual'],differences=[])
        for fraction in (1/8,1/32,1/128,1/512):
            epsilon=original_dt*fraction
            plus=probe.response(state.added(derivative,epsilon),rate=False,time=t)['total_energy']
            minus=probe.response(state.added(derivative,-epsilon),rate=False,time=t)['total_energy']
            row['differences'].append(dict(epsilon_s=epsilon,directional_energy_rate_w=(plus-minus)/(2*epsilon)))
        report['directional_power'].append(row);save()
    for factor in factors:
        _,model,_=build(config);state=initial;steps=(last-first)*factor;dt=(end-start)/steps
        folder=out/f'factor_{factor}';folder.mkdir();energies=[];failure=None
        a=model.response(state,time=start);energy0=a['total_energy']
        energies.append([start,energy0,0.,0.,0.,a['pair_power_residual']])
        for index in range(steps):
            t=start+index*dt
            try:
                previous=model.response(state,time=t)
                if previous['external_power']!=0 or np.any([wave.derivative(t+dt/2)!=0 or wave.derivative(t+dt)!=0 for _,wave in model.source_groups]):
                    raise ValueError('Chosen interval leaves the exact constant-current source hold')
                first_half,_=model.diffusion_step(state,.5*dt)
                e1=model.response(first_half,rate=False,time=t)['total_energy']
                central=model.central_step(first_half,dt,t)
                e2=model.response(central,rate=False,time=t+dt)['total_energy']
                candidate,_=model.diffusion_step(central,.5*dt)
                final=model.response(candidate,time=t+dt);e3=final['total_energy']
            except (ValueError,RuntimeError) as error:
                failure=dict(attempted_step=index+1,type=type(error).__name__,message=str(error));break
            state=candidate
            np.savez_compressed(folder/f'step_{index+1:05}.npz',**arrays(state))
            energies.append([t+dt,e3,e1-previous['total_energy'],e2-e1,e3-e2,final['pair_power_residual']])
        energies=np.asarray(energies);np.savez_compressed(folder/'energy_history.npz',
            time_s=energies[:,0],energy_j=energies[:,1],first_diffusion_change_j=energies[:,2],
            central_change_j=energies[:,3],last_diffusion_change_j=energies[:,4],pair_power_residual_w=energies[:,5])
        entry=dict(factor=factor,dt_s=dt,steps=steps,accepted_steps=len(energies)-1,success=failure is None,failure=failure,
            initial_energy_j=energy0,final_energy_j=float(energies[-1,1]),energy_change_j=float(energies[-1,1]-energy0),
            cumulative_diffusion_energy_change_j=float(energies[:,[2,4]].sum()),
            cumulative_central_energy_change_j=float(energies[:,3].sum()),
            maximum_pair_power_residual_w=float(abs(energies[:,5]).max()))
        if factor==1 and failure is None:
            entry['maximum_original_endpoint_array_differences']={k:float(np.max(abs(v-arrays(saved_final)[k]))) for k,v in arrays(state).items()}
        records.append(entry);endpoints.append(state);save()
        print('Hold interval',factor,entry['success'],entry['energy_change_j'],flush=True)
    if all(row['success'] for row in records):
        fluid,model,meta=build(config);mesh=model.space.mesh;x,w=np.polynomial.legendre.leggauss(10)
        r=(.5*(mesh.rf[1:]+mesh.rf[:-1])[:,None]+.5*np.diff(mesh.rf)[:,None]*x).ravel()
        z=(.5*(mesh.zf[1:]+mesh.zf[:-1])[:,None]+.5*np.diff(mesh.zf)[:,None]*x).ravel()
        rr,zz=np.meshgrid(r,z,indexing='ij')
        dv=2*np.pi*rr*(.5*np.diff(mesh.rf)[:,None]*w).ravel()[:,None]*(.5*np.diff(mesh.zf)[:,None]*w).ravel()[None,:]
        density=meta['central_density_m3']*((1-(rr/mesh.radius)**2)*(1-(zz/mesh.half_length)**2))**1.5
        samples=[material_samples(model,state,rr,zz,density,time=end) for state in endpoints]
        volume=dv*samples[-1]['volume_ratio'];report['endpoint_comparisons_to_finest']=[]
        for factor,sample in zip(factors[:-1],samples[:-1]):
            row=dict(factor=factor,reference_factor=factors[-1])
            for key in ('magnetic','current','velocity'):
                delta=sample[key]-samples[-1][key];reference=samples[-1][key]
                row[key+'_relative_l2']=float(np.sqrt(np.sum(volume[...,None]*delta**2)/np.sum(volume[...,None]*reference**2)))
            row['maximum_position_difference_m']=float(np.linalg.norm(sample['points']-samples[-1]['points'],axis=-1).max())
            report['endpoint_comparisons_to_finest'].append(row)
    fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
    for row in records:
        with np.load(out/f"factor_{row['factor']}"/'energy_history.npz') as data:
            ax.plot((data['time_s']-start)*1e9,(data['energy_j']-row['initial_energy_j'])*1e3,label=f"dt={row['dt_s']*1e9:g} ns")
    ax.axhline(0,color='black',ls='--',label='Independent constant-energy limit')
    ax.set(xlabel='Time from saved hold state [ns]',ylabel='Total energy change [mJ]',title=f"Actual strong-pulse hold: {report['central_integrator']}; zero source work")
    ax.legend();ax.grid(alpha=.2);fig.savefig(out/'hold_energy.png',dpi=170);fig.savefig(out/'hold_energy.pdf');plt.close(fig)
    report['terminal']=True
    if provenance(['examples.material_hold_interval'])!=identity:raise RuntimeError('Interval source changed')
    save();print('Hold continuation terminal',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source');parser.add_argument('destination')
    parser.add_argument('--central-integrator',choices=('rk4','dop853_fixed'))
    parser.add_argument('--factors',type=int,nargs='+',default=[1,2,4])
    args=parser.parse_args();run(args.source,args.destination,central_integrator=args.central_integrator,factors=args.factors)
