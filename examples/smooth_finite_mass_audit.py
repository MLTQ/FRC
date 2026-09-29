"""Read-only independent audits of smooth coupled rate and trajectory records."""
from pathlib import Path
import hashlib
import json
import tarfile
from types import SimpleNamespace
import numpy as np
from plasma_bench.finite_mass_state import evaluate
from plasma_bench.moving_current import MovingCurrent
from .smooth_finite_mass_cases import initial
from .smooth_finite_mass_records import provenance
from .smooth_finite_mass_measurements import contact_rate,contact_changes
from .finite_mass_audit import load,state_from,digest,compare_metrics
from .finite_mass_rate_norms import measure as rate_norms
from .finite_mass_evolution_norms import measure as evolution_norms
from .finite_mass_wave_case import norms as wave_norms
from .finite_mass_contact_case import norms as contact_norms,averages
from .finite_mass_physical_audit import audit as physical
from .finite_mass_manufactured import cell_average,point_conserved
from .moving_current_audit import faraday
from .smooth_induction_norms import measure as induction_norms
from .contact_physical_norms import measure as contact_physical_norms


def source(directory):
    record=json.loads((directory/'source.json').read_text());assert record['provenance']==provenance()
    archive=directory/'source_snapshot.tar.gz';assert digest(archive)==record['archive_sha256']
    with tarfile.open(archive) as tar:
        for name,sha in record['provenance']['source_files'].items():
            assert digest(name)==sha,name
            assert hashlib.sha256(tar.extractfile(name).read()).hexdigest()==sha,name
    return record


def rate(directory,case,frozen):
    record,arrays=load(directory/'rates',case['name']);assert record['configuration']==case and record['provenance']==frozen
    ref,mesh,state,par,time,forcing=initial(case);assert record['parameters']==vars(par)
    for a,b in zip(state.arrays(),state_from(arrays).arrays()):np.testing.assert_array_equal(a,b)
    if not record['completed']:
        assert record['failure'] and 'rate_conserved' not in arrays
        return record,dict(completed=False,record_sha256=digest(directory/'rates'/(case['name']+'.json')))
    ev=evaluate(mesh,state,par);model=ev.model;geo=model.geometry
    result=SimpleNamespace(evaluation=ev,state=state_from(arrays,'rate_'),number_density_rate=arrays['density_rate'])
    if case['family']=='contact':metrics,targets=contact_rate(mesh,ref,state,par,result)
    else:metrics,targets=rate_norms(mesh,ref,state,par,result,time,family=case['family'])
    compare_metrics(metrics,record['metrics'])
    for key,value in targets.items():np.testing.assert_allclose(value,arrays[key],rtol=4e-13,atol=3e-11)
    div=-(np.diff(arrays['radial_flux'],axis=0)+np.diff(arrays['axial_flux'],axis=1))/mesh.volumes[...,None]+arrays['fluid_source']
    for col in (0,2,3,4):np.testing.assert_allclose(div[...,col],result.state.conserved[...,col],rtol=2e-13,atol=1e-16)
    mapping=model.mapping
    nrate=mapping.density(result.state.conserved[...,0])/par.masses.total
    np.testing.assert_allclose(nrate,result.number_density_rate,rtol=2e-14,atol=1e-30)
    p,g=faraday(mesh,arrays['electric']);packed=geo.pack(result.state.psi,result.state.g)
    difference=geo.pack(p,g)-packed;scale=abs(geo.faraday)@abs(arrays['electric'])+abs(packed)
    backward=float(np.linalg.norm(difference)/max(np.linalg.norm(scale),np.finfo(float).tiny))
    assert backward<par.linear_tolerance
    full_relative=float(np.sqrt(difference@(model.hessian@difference))/max(np.sqrt(packed@(model.hessian@packed)),np.finfo(float).tiny))
    terms,potential=MovingCurrent(model).terms(state.psi,state.g,state.electron_entropy,ev.ion_entropy,arrays['regular_velocity'],nrate)
    f=sum(terms.values())+model.mass@model.gradient(potential)
    if forcing is not None:f+=model.mass@geo.project_electric(forcing(time,ev).electric)
    current_rate=model.current(packed);inertia=model.inertia@current_rate
    ohm=model.mass@arrays['electric']-inertia-f
    ohm_backward=float(np.linalg.norm(ohm)/max(np.linalg.norm(abs(model.mass@arrays['electric'])+abs(inertia)+abs(f)),np.finfo(float).tiny))
    assert ohm_backward<par.linear_tolerance
    induction={};uncertainty=0.
    if case['family']!='contact':
        induction=induction_norms(mesh,ref,state,par,arrays,time,family=case['family'])
        higher=induction_norms(mesh,ref,state,par,arrays,time,family=case['family'],order=12)
        uncertainty=max(abs(induction[key]['relative_L2']-higher[key]['relative_L2']) for key in induction)
        assert uncertainty<1e-7
    return record,dict(completed=True,record_sha256=digest(directory/'rates'/(case['name']+'.json')),
        faraday_backward=backward,full_faraday_rate_relative=full_relative,ohm_backward=ohm_backward,
        physical_initial=physical(mesh,state,par,ev),induction=induction,norm_quadrature_absolute_change=uncertainty)


def trajectory(directory,case,frozen):
    record,arrays=load(directory/'trajectories',case['name']);assert record['configuration']==case and record['provenance']==frozen
    ref,mesh,before,par,start,forcing=initial(case);assert record['parameters']==vars(par)
    if record.get('failure_before_trajectory'):
        assert not record['completed'] and record['failure']
        for a,b in zip(before.arrays(),state_from(arrays).arrays()):np.testing.assert_array_equal(a,b)
        return record,dict(completed=False,record_sha256=digest(directory/'trajectories'/(case['name']+'.json')))
    for a,b in zip(before.arrays(),state_from(arrays,'initial_').arrays()):np.testing.assert_array_equal(a,b)
    state=state_from(arrays,'final_');ev=evaluate(mesh,state,par);end=record['actual_time_s'];h=record['history_columns']
    count=len(h['time_s']);assert count>=1 and all(len(v)==count for v in h.values())
    assert h['time_s'][0]==start and h['time_s'][-1]==end and h['steps']==list(range(count))
    assert np.all(np.diff(h['time_s'])>0) and np.all(np.diff(h['time_s'])<=case['dt_s']*(1+2e-12))
    if record['completed']:assert end==record['target_time_s'] and record['failure'] is None
    else:assert record['failure'] and record['attempted_step']['start_s']==end
    for key in ('min_electron_pressure_pa','min_ion_pressure_pa','min_number_density_m3'):assert min(h[key])>0
    budget=np.asarray(h['energy_j'])-h['energy_j'][0]-np.asarray(h['external_energy_input_j'])
    np.testing.assert_allclose(budget,h['energy_budget_residual_j'],rtol=0.,atol=5e-17*abs(h['energy_j'][0]))
    energy_residual=float(np.max(abs(budget))/abs(h['energy_j'][0]));assert energy_residual<2e-13
    drifts={}
    for key,col in (('mass_kg',0),('axial_momentum_kg_m_s',2),('energy_j',3),('angular_momentum_kg_m2_s',4)):
        integral=float(np.sum(mesh.volumes*state.conserved[...,col]));np.testing.assert_allclose(integral,h[key][-1],rtol=2e-15,atol=1e-30)
        if col==0 or (case['family']!='manufactured' and col in (2,4)):
            scale=float(np.sum(mesh.volumes*abs(before.conserved[...,col])))
            # Exact zero angular momentum contacts have zero flux identically;
            # allow only an absolute roundoff scale tied to axial momentum * R.
            if scale==0:scale=float(np.sum(mesh.volumes*abs(before.conserved[...,2])))*mesh.radius
            drift=float(np.max(abs(np.asarray(h[key])-h[key][0])))/max(scale,np.finfo(float).tiny)
            assert drift<2e-13;drifts[key]=drift
    assert max(h['max_linear_residual'])<=par.linear_tolerance
    metrics={};higher={};source_work=None
    if case['family']=='manufactured':
        metrics=evolution_norms(mesh,ref,before,state,par,start,end);higher=evolution_norms(mesh,ref,before,state,par,start,end,order=12)
        energies=[float(np.sum(mesh.volumes*cell_average(mesh,lambda r,z:point_conserved(ref,r,z,t)[...,3],12))) for t in (start,end)]
        difference=energies[1]-energies[0];error=abs(difference-h['external_energy_input_j'][-1])
        assert error<1e-7*max(abs(difference),1e-10)
        source_work=dict(parent_energy_change_j=difference,recorded_work_j=h['external_energy_input_j'][-1],absolute_error_j=error)
    elif case['family']=='wave':metrics=wave_norms(mesh,ref,end,state,par);higher=wave_norms(mesh,ref,end,state,par,order=12)
    elif case['family']=='contact':
        metrics=contact_norms(mesh,ref,state,par,end);higher=contact_norms(mesh,ref,state,par,end,order=18)
        compare_metrics(contact_norms(mesh,ref,before,par,start),record['initial_metrics'])
        assert contact_changes(mesh,ref,before,state,par,start,end)==record['change_metrics']
    compare_metrics(metrics,record['metrics'])
    uncertainty=max((abs(metrics[k]['relative_L2']-higher[k]['relative_L2']) for k in metrics),default=0.);assert uncertainty<1e-7
    initial_ev=evaluate(mesh,before,par)
    a=physical(mesh,before,par,initial_ev);b=physical(mesh,state,par,ev)
    contact_physical={}
    if case['family']=='contact':
        for label,s,e,t in (('initial',before,initial_ev,start),('final',state,ev,end)):
            measured=contact_physical_norms(mesh,ref,s,par,e,t)
            control=contact_physical_norms(mesh,ref,s,par,e,t,order=14)
            change=max(abs(measured['fields'][k]['relative_L2']-control['fields'][k]['relative_L2']) for k in measured['fields'])
            assert change<1e-7
            uncertainty=max(uncertainty,change)
            contact_physical[label]=measured
    return record,dict(completed=record['completed'],record_sha256=digest(directory/'trajectories'/(case['name']+'.json')),
        initial_physical=a,final_physical=b,max_relative_energy_budget=energy_residual,
        normalized_integral_drifts=drifts,source_work=source_work,norm_quadrature_absolute_change=uncertainty,
        contact_physical=contact_physical,
        electron_physical_entropy_change_per_particle=(b['electron_physical_log_entropy_particles']-a['electron_physical_log_entropy_particles'])/a['electron_particles'],
        ion_physical_entropy_change_per_particle=(b['ion_physical_log_entropy_particles']-a['ion_physical_log_entropy_particles'])/a['ion_particles'])
