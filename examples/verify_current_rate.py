"""Audit physical-rate reference convergence, archived fields and measured errors."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plasma_bench.mhd_mesh import MHDMesh
from .electron_reference import ElectronReference
from .current_rate_reference import CurrentRateReference
from .current_inertia_reference import g_rate
from .current_thermal_records import provenance
from .verification_records import write_record


def archive(directory,entry):
    path=directory/entry['fields']['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['fields']['sha256']
    with np.load(path,allow_pickle=False) as z:result={k:z[k] for k in z.files}
    assert all(np.isfinite(v).all() for v in result.values())
    return result


def verify(directory=Path('outputs/current_rate_verification')):
    data=json.loads((directory/'measurements.json').read_text());assert data['completed']
    assert data['provenance']['engine']==provenance()
    for name,sha in data['provenance']['rate_sources'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
    assert [c['nr'] for c in data['cases']]==[8,16,32,64,128]
    assert [c['radial_order'] for c in data['continuum']]==[64,128,256]
    reference=ElectronReference(**data['reference']);reference_audits=[];finest=None
    for entry in data['continuum']:
        a=archive(directory,entry);oracle=CurrentRateReference(reference,entry['radial_order'],entry['axial_points'])
        np.testing.assert_array_equal(a['s'],oracle.s);np.testing.assert_array_equal(a['z'],oracle.z)
        r=reference.radius*np.sqrt(oracle.s[:,None]);r[0]=reference.radius*1e-8
        source={'full':oracle.original_source(),'convective':g_rate(reference,r,oracle.z[None],oracle.mu,oracle.delta)}
        residual={}
        for key in source:
            source[key][-1]=0.
            residual[key]=float(np.linalg.norm(oracle.apply(a[key])-source[key])/np.linalg.norm(source[key]))
            assert residual[key]<2e-11
        reference_audits.append(dict(radial_order=entry['radial_order'],replayed_residual=residual))
        finest=oracle,a
    # Separate axial and solver-tolerance control, compared in physical cell means.
    control=CurrentRateReference(reference,128,128);r=reference.radius*np.sqrt(control.s[:,None]);r[0]=reference.radius*1e-8
    controls={};mesh=MHDMesh(reference.radius,reference.length/2,128,256)
    for key,rhs in (('full',control.original_source()),('convective',g_rate(reference,r,control.z[None],control.mu,control.delta))):
        value,diagnostic=control.solve(rhs,tolerance=1e-12)
        target=finest[0].cell_average(mesh,finest[1][key]);actual=control.cell_average(mesh,value)
        controls[key]=dict(solve=diagnostic,relative_difference=float(np.linalg.norm(actual-target)/np.linalg.norm(target)))
        assert controls[key]['relative_difference']<2e-10
    replay=[];stored={}
    for entry in data['cases']:
        a=archive(directory,entry);stored[entry['nr']]=a
        mesh=MHDMesh(reference.radius,reference.length/2,entry['nr'],entry['nz'])
        np.testing.assert_array_equal(a['radial_faces'],mesh.rf);np.testing.assert_array_equal(a['axial_faces'],mesh.zf)
        assert max(entry['linear_residuals'])<1e-10;values={}
        for key in ('full','convective'):
            target=finest[0].cell_average(mesh,finest[1][key],order=16)
            change=float(np.linalg.norm(target-a[key+'_reference'])/np.linalg.norm(target));assert change<2e-10
            error=a[key+'_computed']-target
            metrics=dict(relative_L2=float(np.sqrt(np.sum(mesh.volumes*error**2)/np.sum(mesh.volumes*target**2))),
                axis_relative_L2=float(np.linalg.norm(error[0])/np.linalg.norm(target[0])),
                maximum_absolute_error=float(abs(error).max()),maximum_reference=float(abs(target).max()))
            for k,v in metrics.items():np.testing.assert_allclose(v,entry['metrics'][key][k],rtol=3e-8,atol=1e-8)
            values[key]=dict(radial_quadrature_relative=change,metrics=metrics)
        replay.append(dict(nr=entry['nr'],values=values))
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
    for key,label in (('full','full physical rate'),('convective','isolated convective contribution')):
        for norm,style in (('relative_L2','o-'),('axis_relative_L2','s--')):
            axes[0].loglog([v['nr'] for v in data['cases']],[v['metrics'][key][norm] for v in data['cases']],style,
                label=label+(', axis' if norm.startswith('axis') else ', volume'))
    axes[0].set(title='After the physical finite-inertia solve',xlabel='radial cells',ylabel='relative gdot error');axes[0].legend(fontsize=8)
    a=stored[128];z=.5*(a['axial_faces'][:-1]+a['axial_faces'][1:])*1e3
    axes[1].plot(z,a['convective_reference'][0],label='independent continuum inversion')
    axes[1].plot(z,a['convective_computed'][0],'--',label='discrete inertia solve')
    axes[1].set(title='First-annulus convective gdot, nr=128',xlabel='z [mm]',ylabel='gdot contribution [T/(m s)]');axes[1].legend(fontsize=8)
    fig.suptitle('Physical magnetic-rate convergence is distinct from intermediate B* forcing',fontsize=13)
    image=directory/'physical_current_rate.png';fig.savefig(image,dpi=160);plt.close(fig)
    names=('examples/verify_current_rate.py','tests/test_current_rate_reference.py')
    report=dict(audit_passed=True,provenance=data['provenance'],reference_audits=reference_audits,
        axial_128_and_tolerance_1e12_controls=controls,replay=replay,
        report_sources={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names},
        figure=dict(file=image.name,sha256=hashlib.sha256(image.read_bytes()).hexdigest()),
        scope='Instantaneous original finite-mass magnetic-rate consistency; intermediate forcing defect remains recorded, and no moving/formation claim is made.')
    write_record(directory/'report.json',report);print(json.dumps(data['cases'][-1]['metrics'],indent=2),flush=True)


if __name__=='__main__':verify()
