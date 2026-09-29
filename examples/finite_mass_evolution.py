"""Reproducible coupled trajectories with complete accepted budgets and failures."""
from dataclasses import asdict,replace
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import numpy as np
from scipy.constants import elementary_charge as E,atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_mass_mhd import FiniteMassMHD
from .moving_species_reference import MovingSpeciesReference
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_manufactured import initialize,ManufacturedSources
from .finite_mass_wave_case import initialize as wave_initial,norms as wave_norms
from .finite_mass_evolution_norms import measure as manufactured_norms
from .finite_mass_records import provenance
from .verification_records import write_record


def protocol():
    ref=MovingSpeciesReference();cases=[]
    def add(name,nr,family,reference,dt=1e-9):
        cases.append(dict(name=name,nr=nr,family=family,reference=asdict(reference),dt_s=dt,
            phase=.23 if family!='wave' else 0.,duration_s=.5/reference.frequency if family!='wave' else .5*np.pi/abs(reference.omega)))
    for nr in (8,16,32):
        add(f'moving_{nr}',nr,'manufactured',ref)
        add(f'target_{nr}',nr,'manufactured',replace(ref,density=1e20))
        add(f'unforced_{nr}',nr,'unforced',ref)
    for label,dt in (('coarse',2e-9),('half',5e-10),('quarter',2.5e-10)):
        add('moving_time_'+label,16,'manufactured',ref,dt)
    wave=FiniteMassWave(.035,.24,1e18,39.948*atomic_mass,.005,6e18*E,gamma_i=1.4,amplitude=1.,branch=2)
    for nr in (8,16,32):add(f'wave_{nr}',nr,'wave',wave)
    add('wave_time_coarse',16,'wave',wave,2e-9);add('wave_time_half',16,'wave',wave,5e-10)
    add('wave_weak_32',32,'wave',replace(wave,amplitude=.1))
    return cases


def case_state(case):
    ref=(FiniteMassWave if case['family']=='wave' else MovingSpeciesReference)(**case['reference'])
    mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],2*case['nr'])
    if case['family']=='wave':state,par=wave_initial(mesh,ref);time=0.;source=None
    else:
        time=case['phase']/ref.frequency;state,par=initialize(mesh,ref,time)
        source=ManufacturedSources(mesh,ref) if case['family']=='manufactured' else None
    return ref,mesh,state,par,time,source


def measure(case):
    ref,mesh,initial,par,start,source=case_state(case);end=start+case['duration_s']
    model=FiniteMassMHD(mesh,initial,par,max_timestep=case['dt_s'],source=source,initial_time=start)
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
    arrays={}
    for prefix,state in (('initial',initial),('final',model.state)):
        for name,value in zip(('conserved','psi','g','electron_entropy'),state.arrays()):arrays[prefix+'_'+name]=value
    metrics={}
    if case['family']=='manufactured':metrics=manufactured_norms(mesh,ref,initial,model.state,par,start,model.time)
    elif case['family']=='wave':metrics=wave_norms(mesh,ref,model.time,model.state,par)
    record=dict(configuration=case,parameters=asdict(par),completed=failure is None and model.time==end,failure=failure,
        attempted_step=attempted if failure else None,target_time_s=end,actual_time_s=model.time,
        metrics=metrics,history_columns={key:[row[key] for row in history] for key in history[0]},
        scope='Complete coupled RK4 trajectory in a positive-density periodic cylinder; no physical exterior, transport or formation validation.')
    return record,arrays


def run(directory,group=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);frozen=provenance()
    if group is not None and group not in {c['name'] for c in protocol()}:raise ValueError('Unknown case')
    for case in protocol():
        if group is not None and case['name']!=group:continue
        start=perf_counter();record,arrays=measure(case);assert provenance()==frozen
        path=directory/(case['name']+'.npz');np.savez_compressed(path,**arrays)
        record.update(provenance=frozen,elapsed_s=perf_counter()-start,fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        write_record(directory/(case['name']+'.json'),record)
        print(case['name'],'completed',record['completed'],'time',record['actual_time_s'],'failure',record['failure'],'metrics',record['metrics'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/finite_mass_evolution');parser.add_argument('--group')
    args=parser.parse_args();run(args.directory,args.group)
