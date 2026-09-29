"""Read-only audits of the frozen complete-rate and actual evolution records."""
from pathlib import Path
import hashlib
import json
import tarfile
from types import SimpleNamespace
import numpy as np
from plasma_bench.finite_mass_state import FiniteMassState,FiniteMassParameters,evaluate
from .finite_mass_rate_measurements import case_state as rate_initial,protocol as rate_protocol
from .finite_mass_rate_norms import measure as rate_norms
from .finite_mass_evolution import case_state as evolution_initial
from .finite_mass_evolution_norms import measure as trajectory_norms
from .finite_mass_wave_case import norms as wave_norms
from .finite_mass_physical_audit import audit as physical
from .finite_mass_manufactured import cell_average,point_conserved
from .moving_current_audit import faraday


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(directory,name):
    path=directory/(name+'.json');record=json.loads(path.read_text());field=directory/record['fields']['file']
    assert field.name==name+'.npz' and digest(field)==record['fields']['sha256']
    with np.load(field,allow_pickle=False) as file:arrays={k:file[k] for k in file.files}
    assert all(np.isfinite(a).all() for a in arrays.values())
    return record,arrays


def state_from(arrays,prefix=''):
    return FiniteMassState(*(arrays[prefix+key] for key in ('conserved','psi','g','electron_entropy')))


def compare_metrics(actual,stored,rtol=5e-12):
    assert set(actual)==set(stored)
    for key,metrics in actual.items():
        for name,value in metrics.items():np.testing.assert_allclose(value,stored[key][name],rtol=rtol,atol=2e-12)


def baseline(directory):
    manifest=json.loads((directory/'baseline.json').read_text());source=manifest['provenance']
    assert digest(directory/'source_snapshot.tar.gz')==manifest['source_archive_sha256']
    with tarfile.open(directory/'source_snapshot.tar.gz') as tar:
        for path,sha in source['source_files'].items():assert hashlib.sha256(tar.extractfile(path).read()).hexdigest()==sha
    # Only the provenance driver's expanded roots differ since this rate snapshot.
    for path,sha in source['source_files'].items():
        if path!='examples/finite_mass_records.py':assert digest(path)==sha,path
    records=[];audits={}
    for case in rate_protocol():
        record,arrays=load(directory,case['name']);assert record['configuration']==case and record['provenance']==source
        ref,mesh,initial,par,forcing,time=rate_initial(case);state=state_from(arrays)
        assert record['parameters']==vars(par) and record['input_unchanged']
        for a,b in zip(state.arrays(),initial.arrays()):np.testing.assert_array_equal(a,b)
        assert record['completed'],record['failure']
        ev=evaluate(mesh,state,par);rate=SimpleNamespace(evaluation=ev,state=state_from(arrays,'rate_'),number_density_rate=arrays['density_rate'])
        metrics,targets=rate_norms(mesh,ref,state,par,rate,time,family=case['family'])
        compare_metrics(metrics,record['metrics'])
        for key,value in targets.items():np.testing.assert_allclose(value,arrays[key],rtol=3e-13,atol=3e-11)
        refined,_=rate_norms(mesh,ref,state,par,rate,time,family=case['family'],order=12)
        uncertainty=max(abs(refined[k]['relative_L2']-metrics[k]['relative_L2']) for k in metrics)
        assert uncertainty<1e-7
        volume=mesh.volumes[...,None];div=-(np.diff(arrays['radial_flux'],axis=0)+np.diff(arrays['axial_flux'],axis=1))/volume+arrays['fluid_source']
        for column in (0,2,3,4):np.testing.assert_allclose(div[...,column],rate.state.conserved[...,column],rtol=2e-13,atol=1e-16)
        fp,fg=faraday(mesh,arrays['electric']);packed=ev.model.geometry.pack(fp-rate.state.psi,fg-rate.state.g)
        exact=ev.model.geometry.pack(rate.state.psi,rate.state.g);h=ev.model.hessian
        ferr=float(np.sqrt(packed@(h@packed)/(exact@(h@exact))))
        assert ferr<5*par.linear_tolerance
        audits[case['name']]=dict(record_sha256=digest(directory/(case['name']+'.json')),norm_quadrature_absolute_change=uncertainty,faraday_relative=ferr,
            physical_initial=physical(mesh,state,par,ev))
        records.append(record);print('rate audit',case['name'],flush=True)
    return manifest,records,audits


def evolution(directory,case,source):
    record,arrays=load(directory,case['name']);assert record['configuration']==case and record['provenance']==source
    ref,mesh,initial,par,start,forcing=evolution_initial(case);assert record['parameters']==vars(par)
    for a,b in zip(initial.arrays(),state_from(arrays,'initial_').arrays()):np.testing.assert_array_equal(a,b)
    state=state_from(arrays,'final_');ev=evaluate(mesh,state,par);end=record['actual_time_s']
    h=record['history_columns'];count=len(h['time_s']);assert all(len(v)==count for v in h.values())
    assert count>=1 and h['time_s'][0]==start and h['time_s'][-1]==end
    assert np.all(np.diff(h['time_s'])>0) and np.all(np.diff(h['time_s'])<=case['dt_s']*(1+1e-12))
    assert h['steps']==list(range(count))
    for key in ('min_electron_pressure_pa','min_ion_pressure_pa','min_number_density_m3'):assert min(h[key])>0
    if record['completed']:assert record['failure'] is None and end==record['target_time_s']
    else:assert record['failure'] and record['attempted_step']['start_s']==end
    energies=np.array(h['energy_j']);work=np.array(h['external_energy_input_j'])
    budget=energies-energies[0]-work
    np.testing.assert_allclose(budget,h['energy_budget_residual_j'],rtol=0.,atol=5e-17*abs(energies[0]))
    assert np.max(abs(budget))/abs(energies[0])<2e-13
    for key,column in (('mass_kg',0),('axial_momentum_kg_m_s',2),('energy_j',3),('angular_momentum_kg_m2_s',4)):
        value=float(np.sum(mesh.volumes*state.conserved[...,column]));np.testing.assert_allclose(value,h[key][-1],rtol=2e-15,atol=1e-30)
    drifts={}
    for key,column in (('mass_kg',0),('axial_momentum_kg_m_s',2),('angular_momentum_kg_m2_s',4)):
        if key!='mass_kg' and case['family']=='manufactured':continue
        # Oscillatory net momentum can vanish: normalize by the initial integral
        # of its absolute local density, rather than divide by a cancelled sum.
        scale=float(np.sum(mesh.volumes*abs(initial.conserved[...,column])))
        drift=float(np.max(abs(np.asarray(h[key])-h[key][0])))
        relative=drift/max(scale,np.finfo(float).tiny)
        assert relative<2e-13
        drifts[key]=dict(maximum_absolute_change=drift,initial_absolute_integral=scale,normalized_drift=relative)
    assert max(h['max_linear_residual'])<=par.linear_tolerance
    metrics={};refined={};source_energy=None
    if case['family']=='manufactured':
        metrics=trajectory_norms(mesh,ref,initial,state,par,start,end)
        refined=trajectory_norms(mesh,ref,initial,state,par,start,end,order=12)
        continuum=[]
        for time in (start,end):continuum.append(float(np.sum(mesh.volumes*cell_average(mesh,lambda r,z:point_conserved(ref,r,z,time)[...,3],12))))
        source_energy=dict(exact_parent_energy_change_j=continuum[1]-continuum[0],recorded_source_work_j=work[-1],absolute_difference_j=abs(work[-1]-(continuum[1]-continuum[0])))
        # The declared RK quadrature has finite time error, retained separately.
        assert source_energy['absolute_difference_j']<1e-7*max(abs(continuum[1]-continuum[0]),1e-10)
    elif case['family']=='wave':
        metrics=wave_norms(mesh,ref,end,state,par);refined=wave_norms(mesh,ref,end,state,par,order=12)
    compare_metrics(metrics,record['metrics'])
    uncertainty=max((abs(metrics[k]['relative_L2']-refined[k]['relative_L2']) for k in metrics),default=0.)
    assert uncertainty<1e-7
    a=physical(mesh,initial,par,evaluate(mesh,initial,par));b=physical(mesh,state,par,ev)
    return record,dict(record_sha256=digest(directory/(case['name']+'.json')),initial_physical=a,final_physical=b,
        max_relative_energy_budget=float(np.max(abs(budget))/abs(energies[0])),source_energy=source_energy,
        norm_quadrature_absolute_change=uncertainty,unforced_integral_drifts=drifts,
        electron_physical_entropy_change_per_particle=(b['electron_physical_log_entropy_particles']-a['electron_physical_log_entropy_particles'])/a['electron_particles'],
        ion_physical_entropy_change_per_particle=(b['ion_physical_log_entropy_particles']-a['ion_physical_log_entropy_particles'])/a['ion_particles'])
