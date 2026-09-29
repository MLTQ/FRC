"""Strict archival and independent physical energy audit; no state modifications."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import electron_mass as ME,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from .current_thermal_cases import protocol
from .current_thermal_records import provenance
from .current_beltrami_reference import CurrentBeltrami
from .electron_reference import ElectronReference
from .current_thermal_norms import physical_energy,wave_norms,wave_face_power


def endpoint(mesh,case,archive,prefix,diagnostic):
    n=archive['number_density'];ue=archive[prefix+'_electron_entropy'];ui=archive[prefix+'_ion_entropy']
    psi,g,current=(archive[prefix+'_'+k] for k in ('psi','g','current'))
    me=ME if case['model']=='two_fluid' else 0.;mu=me*case['ion_mass']/(me+case['ion_mass'])
    energies=[]
    for order in (8,12):
        energies.append(physical_energy(mesh,psi,g,current,n,ue,ui,reduced_mass=mu,
            density_scale=case['density_scale'],gamma_e=case['gamma_e'],gamma_i=case['gamma_i'],order=order))
    replay={};quad={};totals={}
    for key,value in energies[1].items():
        total=float(np.sum(mesh.volumes*value));totals[key+'_energy_j']=total
        denominator=max(np.linalg.norm(value),np.finfo(float).tiny)
        quad[key]=float(np.linalg.norm(value-energies[0][key])/denominator)
        target=diagnostic[key+'_energy_j']
        replay[key]=abs(total-target)/max(abs(total),np.finfo(float).tiny)
        # Positive physical rational quadrature (current/n) differs from order 5;
        # this bound checks the declared smooth protocol, not arbitrary densities.
        assert replay[key]<2e-8,(case['name'],prefix,key,replay[key])
        assert quad[key]<1e-11,(case['name'],prefix,key,quad[key])
    total=archive[prefix+'_total_energy'];ion=total-energies[1]['magnetic']-energies[1]['current_kinetic']-energies[1]['electron']
    defect=(case['gamma_i']-1)*(total-sum(energies[1].values()))
    # Positive quadrature accumulation and endpoint energy subtraction lose
    # absolute digits under a large magnetic background. Record that uncertainty.
    operations=2*(8**2+12**2+case['quadrature']**2)+64
    gamma=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    pressure_roundoff=float(gamma*(case['gamma_i']-1)*np.max(abs(total)+sum(abs(v) for v in energies[1].values())))
    for key,value in dict(conservative_energy_j=float(np.sum(mesh.volumes*total)),
        min_inferred_ion_pressure_pa=float((case['gamma_i']-1)*ion.min()),
        local_pressure_defect_rms_pa=float(np.sqrt(np.sum(mesh.volumes*defect**2)/np.sum(mesh.volumes))),
        local_pressure_defect_max_pa=float(abs(defect).max()),axis_pressure_defect_max_pa=float(abs(defect[0]).max())).items():
        np.testing.assert_allclose(value,diagnostic[key],rtol=2e-7,
            atol=1e-10+(pressure_roundoff if key.endswith('_pa') else 0.),err_msg=key)
    # Dual volumes computed independently as quarter-cell supports.
    volume=np.zeros_like(n)
    for di in (0,1):
        volume[di:di+mesh.nr]+=.5*mesh.volumes
    np.testing.assert_allclose(volume,archive['dual_volumes'],rtol=8*mesh.nr*np.finfo(float).eps)
    for label,u in (('electron',ue),('ion',ui)):
        gamma=case['gamma_e' if label=='electron' else 'gamma_i']
        pressure=(n/case['density_scale'])**(gamma-1)*u
        assert np.all(u>0)
        for key,value in {label+'_entropy_proxy_j':np.sum(volume*u),
                          label+'_log_entropy_particles':np.sum(volume*n*np.log(u/(n*E))),
                          'min_'+label+'_pressure_pa':pressure.min()}.items():
            np.testing.assert_allclose(value,diagnostic[key],rtol=3e-12,atol=1e-20)
    np.testing.assert_array_equal(psi[:,0],psi[:,-1]);np.testing.assert_allclose(psi[0],0,atol=1e-25)
    np.testing.assert_allclose(psi[-1],psi[-1,0],atol=1e-20)
    return dict(energy_replay_relative=replay,quadrature_relative=quad,physical_energy_j=totals,
                pressure_replay_roundoff_bound_pa=pressure_roundoff)


def audit_case(case,record,a,mesh):
    history=record['history_columns'];count=record['accepted_steps'];time=np.asarray(history['time_s'])
    assert all(len(v)==count+1 and np.isfinite(v).all() for v in history.values())
    assert time[0]==0 and np.all(np.diff(time)>0)
    assert all(history[k][0]==v for k,v in record['initial'].items())
    assert all(history[k][-1]==v for k,v in record['final'].items())
    end=record['requested_end_time_s'];dt=end/case['steps']
    np.testing.assert_allclose(time,np.arange(count+1)*dt,rtol=2e-15,atol=0.)
    assert min(history['min_inferred_ion_pressure_pa'])>0
    assert max(history['nonlinear_residual'])<=case['tolerance']
    assert max(history['linear_residual'])<=max(case['tolerance'],1e-12)
    if 'stored_state_residual' in history:
        assert max(history['stored_state_residual'])<=case['tolerance']
    if record['completed']:
        assert count==case['steps'] and abs(time[-1]-end)<=2*np.finfo(float).eps*end and record['failure'] is None
    else:
        assert count<case['steps'] and record['failure']['input_unchanged']
        assert record['failure']['trial_time_s']==(count+1)*dt
    result={prefix:endpoint(mesh,case,a,prefix,record[prefix]) for prefix in ('initial','final')}
    if 'rejected_psi' in a:
        result['rejected']=endpoint(mesh,case,a,'rejected',record['failure']['diagnostic'])
        assert record['failure']['diagnostic']['min_inferred_ion_pressure_pa']<=0
    for key in ('total_energy_j','conservative_energy_j','electron_entropy_proxy_j','ion_entropy_proxy_j'):
        values=np.asarray(history[key]);result['max_'+key+'_relative_drift']=float(np.max(abs(values-values[0]))/abs(values[0]))
    particles=float(np.sum(a['dual_volumes']*a['number_density']))
    for species in ('electron','ion'):
        values=np.asarray(history[species+'_log_entropy_particles'])
        result[species+'_max_log_entropy_change_per_particle']=float(np.max(abs(values-values[0]))/particles)
    if count:
        fr,fz=a['last_radial_power'],a['last_axial_power']
        np.testing.assert_array_equal(fr[0],0.);np.testing.assert_allclose(fr[-1],0.,atol=1e-13)
        np.testing.assert_array_equal(fz[:,0],fz[:,-1])
        # Independently use physical coefficient Faraday, including curl-free E.
        electric=a['last_electric'];cells=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
        er=electric[:cells].reshape(mesh.nr,mesh.nz)
        et=np.pad(electric[cells:cells+theta].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
        ez=electric[cells+theta:].reshape(mesh.nr+1,mesh.nz)
        dg=np.diff(ez,axis=0)/mesh.radial_measure[:,None]-(np.roll(er,-1,axis=1)-er)/mesh.dz
        dp=-mesh.rf[:,None]*et
        discrepancies=[];rounding=[]
        for rate,final,mid in ((dg,a['final_g'],a['last_mid_g']),
                               (dp,a['final_psi'][:,:-1],a['last_mid_psi'][:,:-1])):
            actual=2*(final-mid)/dt
            denominator=max(np.linalg.norm(actual),1e-15)
            absolute=float(np.linalg.norm(rate-actual));relative=absolute/denominator
            # Rounding of the state addition, midpoint, subtraction and division.
            # This bounds recoverable rates; it adds no physical-rate floor.
            bound=8*np.finfo(float).eps*np.linalg.norm((abs(final)+abs(mid))/dt)
            discrepancies.append(relative);rounding.append(dict(absolute=absolute,state_roundoff_bound=float(bound)))
            assert relative<2e-7 or absolute<=bound,(case['name'],relative,bound,absolute)
        result['last_full_electric_faraday_relative']=discrepancies
        result['last_full_electric_faraday_roundoff']=rounding
    if case['kind']=='wave':
        ref=CurrentBeltrami(**case['reference'])
        psi,g=ref.initial_fields(mesh,order=24)
        np.testing.assert_allclose(a['initial_psi'],psi,rtol=1e-14,atol=1e-20)
        np.testing.assert_allclose(a['initial_g'],g,rtol=2e-12,atol=1e-16)
        metrics=[wave_norms(mesh,ref,time[-1],a['final_psi'],a['final_g'],a['final_current'],order=o) for o in (8,12)]
        result['wave_quadrature_absolute']={k:abs(metrics[0][k]-metrics[1][k]) for k in metrics[0]}
        for key,value in metrics[1].items():np.testing.assert_allclose(value,record['metrics'][key],rtol=2e-9,atol=1e-15)
        if count:
            powers=[wave_face_power(mesh,ref,time[-1]-.5*dt,order=o) for o in (10,16)]
            result['face_quadrature_relative']=[float(np.linalg.norm(p-q)/np.linalg.norm(q)) for p,q in zip(*powers)]
    else:
        ref=ElectronReference(**case['reference']);f=ref.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))
        np.testing.assert_array_equal(a['number_density'],f['density'])
        np.testing.assert_allclose((a['number_density']/case['density_scale'])**(case['gamma_e']-1)*a['initial_electron_entropy'],f['pressure'],rtol=3e-16)
        np.testing.assert_allclose(a['initial_g'],ref.cell_average(mesh,'g',order=16),rtol=1e-12,atol=1e-18)
    return result


def collection(directory):
    directory=Path(directory);data=json.loads((directory/'measurements.json').read_text())
    assert data['protocol_finished'] and data['provenance']==provenance()
    assert [r['configuration'] for r in data['cases']]==protocol()
    records={};archives={};meshes={};audits={}
    for case,record in zip(protocol(),data['cases']):
        name=case['name'];assert record==json.loads((directory/(name+'.json')).read_text())
        assert record['provenance']==data['provenance'];path=directory/record['fields']['file']
        assert path.name==name+'.npz' and hashlib.sha256(path.read_bytes()).hexdigest()==record['fields']['sha256']
        with np.load(path,allow_pickle=False) as file:a={key:file[key] for key in file.files}
        assert all(np.isfinite(v).all() for v in a.values())
        mesh=MHDMesh(case['reference']['radius'],case['reference']['length']/2,case['nr'],case['nz'])
        for key,value in (('radial_faces',mesh.rf),('axial_faces',mesh.zf),('primal_volumes',mesh.volumes)):
            np.testing.assert_array_equal(a[key],value)
        audits[name]=audit_case(case,record,a,mesh);records[name]=record;archives[name]=a;meshes[name]=mesh
    return records,archives,meshes,audits
