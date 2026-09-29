"""Read-only archive, independent Faraday and physical-norm moving-rate audit."""
import hashlib
import json
import tarfile
from types import SimpleNamespace
import numpy as np
from .moving_current_case import initialize,norms
from .moving_current_work import verify as verify_work
from plasma_bench.moving_current import MovingCurrent


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def faraday(mesh,electric):
    """Direct physical -curl E coefficients; no production curl/gradient maps."""
    cells=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    er=electric[:cells].reshape(mesh.nr,mesh.nz)
    et=electric[cells:cells+theta].reshape(mesh.nr-1,mesh.nz)
    ez=electric[cells+theta:].reshape(mesh.nr+1,mesh.nz)
    psi=np.zeros((mesh.nr+1,mesh.nz+1));psi[1:-1,:-1]=-mesh.rf[1:-1,None]*et
    psi[:,-1]=psi[:,0]
    g=np.diff(ez,axis=0)/mesh.radial_measure[:,None]-(np.roll(er,-1,axis=1)-er)/mesh.dz
    return psi,g


def audit_case(directory,case,source):
    path=directory/(case['name']+'.json');record=json.loads(path.read_text())
    assert record['configuration']==case and record['provenance']==source
    assert record['fields']['file']==case['name']+'.npz'
    assert digest(directory/record['fields']['file'])==record['fields']['sha256']
    with np.load(directory/record['fields']['file'],allow_pickle=False) as file:
        arrays={k:file[k] for k in file.files}
    assert all(np.isfinite(a).all() for a in arrays.values())
    ref,model,psi,g,ue,ui,v,nt,forcing,time=initialize(case);geo=model.geometry;mesh=model.mesh
    expected=dict(psi=psi,g=g,number_density=model.electron.density,electron_entropy=ue,ion_entropy=ui,
        regular_velocity=v,density_rate=nt,radial_faces=mesh.rf,axial_faces=mesh.zf)
    for key,value in expected.items():np.testing.assert_array_equal(arrays[key],value)
    assert record['input_unchanged']
    if not record['completed']:
        assert record['status']=='rejected' and record['failure'] and set(arrays)==set(expected)
        return record,dict(initial_fields_checked=True,rejection_retained=True,record_sha256=digest(path))
    assert record['status']=='succeeded' and record['failure'] is None
    assert record['linear_residual']<=case['tolerance']
    assert max(record['residual_components'].values())==record['linear_residual']
    rate=geo.pack(arrays['rate_psi'],arrays['rate_g'])
    norm=lambda value:float(np.sqrt(max(value@(model.hessian@value),0.)))
    computed=geo.pack(*faraday(mesh,arrays['electric_full']))
    faraday_error=norm(computed-rate)/norm(rate)
    assert faraday_error<5*case['tolerance'],faraday_error
    np.testing.assert_array_equal(arrays['rate_psi'][[0,-1]],0.)
    np.testing.assert_array_equal(arrays['rate_psi'][:,0],arrays['rate_psi'][:,-1])
    current=model.current(rate)
    np.testing.assert_allclose(arrays['rate_current'],current,rtol=2e-12,atol=2e-12*np.max(abs(current)))
    weak,potential=MovingCurrent(model).terms(psi,g,ue,ui,v,nt)
    np.testing.assert_allclose(arrays['potential'],potential,rtol=2e-14,atol=1e-15)
    full=arrays['electric_reduced']+model.gradient(potential)
    np.testing.assert_allclose(arrays['electric_full'],full,rtol=2e-14,atol=1e-12)
    forcing=sum(weak.values())+model.mass@geo.project_electric(forcing)
    work=verify_work(model,psi,g,ue,ui,nt,SimpleNamespace(magnetic=rate,current=current),forcing)
    for key,value in work.items():
        if isinstance(value,list):
            for row,stored in zip(value,record['work'][key]):
                for field in row:np.testing.assert_allclose(row[field],stored[field],rtol=2e-4,atol=2e-11)
        else:np.testing.assert_allclose(value,record['work'][key],rtol=2e-9,atol=2e-11)
    metrics=norms(mesh,ref,time,arrays['rate_psi'],arrays['rate_g'],arrays['rate_current'],arrays['electric_full'])
    refined=norms(mesh,ref,time,arrays['rate_psi'],arrays['rate_g'],arrays['rate_current'],arrays['electric_full'],order=12)
    uncertainty={}
    for key in metrics:
        for field in metrics[key]:np.testing.assert_allclose(metrics[key][field],record['metrics'][key][field],rtol=3e-13,atol=1e-14)
        uncertainty[key]={field:abs(refined[key][field]-metrics[key][field])
                          for field in ('relative_L2','axis_relative_L2')}
        assert max(uncertainty[key].values())<1e-7,uncertainty
    omitted={}
    for name,stored in record['omission_controls'].items():
        alternative=arrays['without_'+name+'_rate'];change=norm(alternative-rate)/norm(rate)
        np.testing.assert_allclose(change,stored['magnetic_rate_relative_change'],rtol=2e-13)
        jdot=model.current(alternative);rhs=forcing-weak[name]+model.inertia@jdot
        electric=np.concatenate([solve(part) for solve,part in zip(geo.solves,geo.split_electric(rhs))])
        error=norm(geo.pack(*faraday(mesh,electric))-alternative)/norm(alternative)
        assert error<5*case['tolerance']
        omitted[name]=dict(magnetic_relative_change=change,replayed_faraday_relative=error)
    return record,dict(initial_fields_checked=True,record_sha256=digest(path),
        full_electric_faraday_relative=faraday_error,norm_quadrature_absolute_change=uncertainty,
        order12_metrics=refined,energy_work=work,omission_controls=omitted)


def audit_initial(directory):
    """Verify the failed source snapshot without relabeling its results as current."""
    record=json.loads((directory/'attempt.json').read_text());archive=directory/'source_snapshot.tar.gz'
    assert not record['completed'] and digest(archive)==record['source_archive_sha256']
    with tarfile.open(archive) as tar:
        for path,sha in record['provenance']['source_files'].items():
            assert hashlib.sha256(tar.extractfile(path).read()).hexdigest()==sha
    cases=[]
    for path in sorted(directory.glob('*.json')):
        if path.name=='attempt.json':continue
        case=json.loads(path.read_text());assert case['provenance']==record['provenance']
        assert digest(directory/case['fields']['file'])==case['fields']['sha256']
        cases.append(dict(name=path.stem,completed=case['completed'],record_sha256=digest(path)))
    failed=json.loads((directory/'moving_128.json').read_text())
    assert not failed['completed'] and 'residual' in failed['failure']
    assert 'AssertionError' in (directory/'main.log').read_text()
    return dict(attempt=record,cases=cases,logs={p.name:digest(p) for p in directory.glob('*.log')})
