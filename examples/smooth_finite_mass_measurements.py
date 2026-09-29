"""Complete smooth coupled rates/trajectories, including rejected-step evidence."""
from pathlib import Path
import argparse
import hashlib
import json
import tarfile
from time import perf_counter
import numpy as np
from plasma_bench.smooth_finite_mass_mhd import SmoothFiniteMassMHD
from plasma_bench.smooth_finite_mass_rhs import rhs
from plasma_bench.finite_mass_rhs import pressure_rates
from plasma_bench.finite_mass_state import evaluate
from .smooth_finite_mass_cases import initial,rate_protocol,evolution_protocol
from .smooth_finite_mass_records import provenance
from .finite_mass_rate_norms import measure as rate_norms
from .finite_mass_evolution_norms import measure as evolution_norms
from .finite_mass_wave_case import norms as wave_norms
from .finite_mass_contact_case import averages,norms as contact_norms
from .verification_records import write_record


def freeze(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);source=provenance()
    manifest=directory/'source.json'
    if manifest.exists():assert json.loads(manifest.read_text())['provenance']==source
    else:
        path=directory/'source_snapshot.tar.gz'
        with tarfile.open(path,'w:gz') as tar:
            for name in source['source_files']:tar.add(name,arcname=name)
        write_record(manifest,dict(provenance=source,archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return source


def fields(state,prefix=''):
    return {prefix+k:v.copy() for k,v in zip(('conserved','psi','g','electron_entropy'),state.arrays())}


def contact_rate(mesh,ref,state,par,rate):
    exact=averages(mesh,ref,rate=True);pe,pi=pressure_rates(state,rate,par)
    pe_ref=np.broadcast_to(ref.fields(mesh.zf[:-1],rate=True)['electron_pressure'][None],pe.shape)
    items=dict(mass=(rate.state.conserved[...,0],exact['density'],mesh.volumes),
        axial_momentum=(rate.state.conserved[...,2],exact['momentum'],mesh.volumes),
        total_energy=(rate.state.conserved[...,3],exact['energy'],mesh.volumes),
        electron_pressure=(pe,pe_ref,rate.evaluation.model.mapping.volumes),
        ion_pressure=(pi,exact['ion_pressure'],mesh.volumes))
    metrics={key:dict(relative_L2=float(np.sqrt(np.sum(v*(a-b)**2)/np.sum(v*b*b))),
        rms_error=float(np.sqrt(np.sum(v*(a-b)**2)/v.sum())),maximum_error=float(abs(a-b).max())) for key,(a,b,v) in items.items()}
    targets={key+'_target':b for key,(a,b,v) in items.items()}
    targets.update(electron_pressure_rate=pe,ion_pressure_rate=pi)
    return metrics,targets


def contact_changes(mesh,ref,before,after,par,start,end):
    a,b=evaluate(mesh,before,par),evaluate(mesh,after,par)
    nodes=mesh.zf[:-1];exact0,exact1=averages(mesh,ref,start),averages(mesh,ref,end)
    values=dict(electron_pressure=(a.model.electron.pressure(before.electron_entropy),b.model.electron.pressure(after.electron_entropy),
        np.broadcast_to(ref.fields(nodes,start)['electron_pressure'][None],before.electron_entropy.shape),
        np.broadcast_to(ref.fields(nodes,end)['electron_pressure'][None],after.electron_entropy.shape),a.model.mapping.volumes),
        ion_pressure=((par.gamma_i-1)*a.components['ion'],(par.gamma_i-1)*b.components['ion'],exact0['ion_pressure'],exact1['ion_pressure'],mesh.volumes))
    result={}
    for key,(x,y,ea,eb,v) in values.items():
        norm=lambda value:float(np.sqrt(np.sum(v*value*value)))
        physical=norm(eb-ea);error=norm((y-x)-(eb-ea))
        result[key]=dict(error_in_change_rms=error/np.sqrt(v.sum()),physical_change_rms=physical/np.sqrt(v.sum()),
            error_in_change_relative_L2=error/physical if physical>1e-12*norm(ea) else None)
    return result


def rate(case):
    ref,mesh,state,par,time,source=initial(case);saved=fields(state);arrays=fields(state)
    result=rhs(mesh,state,par,time=time,source=source,entropy_transport=case['entropy_transport'])
    if case['family']=='contact':metrics,targets=contact_rate(mesh,ref,state,par,result)
    else:metrics,targets=rate_norms(mesh,ref,state,par,result,time,family=case['family'])
    arrays.update(fields(result.state,'rate_'),electric=result.electric,radial_flux=result.radial_flux,
        axial_flux=result.axial_flux,fluid_source=result.fluid_source,density_rate=result.number_density_rate,
        regular_velocity=result.regular_velocity,**targets)
    for key,value in saved.items():np.testing.assert_array_equal(value,fields(state)[key])
    return dict(completed=True,failure=None,parameters=vars(par),metrics=metrics,
        linear_residual=result.linear_residual,linear_iterations=result.linear_iterations,
        linear_refinements=result.linear_refinements,induction_components=result.induction_components),arrays


def trajectory(case):
    ref,mesh,state,par,start,source=initial(case);end=start+case['duration_s']
    model=SmoothFiniteMassMHD(mesh,state,par,max_timestep=case['dt_s'],source=source,initial_time=start,
        entropy_transport=case['entropy_transport'])
    history=[model.diagnostics()];failure=None;attempted=None
    while model.time<end:
        remaining=end-model.time
        if remaining<=4*np.finfo(float).eps*abs(end):model.time=end;history[-1]['time_s']=end;break
        dt=min(remaining,case['dt_s']);attempted=dict(start_s=model.time,dt_s=dt)
        saved=model.state.copy();before=model.diagnostics()
        try:model.step(dt)
        except (RuntimeError,ValueError) as error:
            failure=repr(error)
            for a,b in zip(saved.arrays(),model.state.arrays()):np.testing.assert_array_equal(a,b)
            assert model.diagnostics()==before
            break
        history.append(model.diagnostics())
    metrics={};initial_metrics={};changes={}
    if case['family']=='manufactured':metrics=evolution_norms(mesh,ref,state,model.state,par,start,model.time)
    elif case['family']=='wave':metrics=wave_norms(mesh,ref,model.time,model.state,par)
    elif case['family']=='contact':
        metrics=contact_norms(mesh,ref,model.state,par,model.time);initial_metrics=contact_norms(mesh,ref,state,par,start)
        changes=contact_changes(mesh,ref,state,model.state,par,start,model.time)
    record=dict(parameters=vars(par),completed=bool(failure is None and model.time==end),failure=failure,
        attempted_step=attempted if failure else None,target_time_s=end,actual_time_s=model.time,
        metrics=metrics,initial_metrics=initial_metrics,change_metrics=changes,
        history_columns={key:[row[key] for row in history] for key in history[0]})
    return record,dict(**fields(state,'initial_'),**fields(model.state,'final_'))


def run(directory,kind,*,partition=0,partitions=1):
    directory=Path(directory);frozen=freeze(directory);target=directory/kind;target.mkdir(exist_ok=True)
    cases=rate_protocol() if kind=='rates' else evolution_protocol()
    for index,case in enumerate(cases):
        if index%partitions!=partition:continue
        begin=perf_counter()
        try:record,arrays=(rate if kind=='rates' else trajectory)(case)
        except (RuntimeError,ValueError) as error:
            # Initialization/rate failure retains the exact positive input when available.
            ref,mesh,state,par,time,source=initial(case)
            record=dict(completed=False,failure=repr(error),parameters=vars(par),metrics={},failure_before_trajectory=True)
            arrays=fields(state)
        assert provenance()==frozen
        path=target/(case['name']+'.npz');np.savez_compressed(path,**arrays)
        record.update(configuration=case,provenance=frozen,elapsed_s=perf_counter()-begin,
            fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        write_record(target/(case['name']+'.json'),record)
        print(kind,case['name'],record['completed'],record['failure'],record.get('actual_time_s'),record['metrics'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/smooth_finite_mass_verification')
    parser.add_argument('--kind',choices=['rates','trajectories'],required=True)
    parser.add_argument('--partition',type=int,default=0);parser.add_argument('--partitions',type=int,default=1)
    args=parser.parse_args();run(args.directory,args.kind,partition=args.partition,partitions=args.partitions)
