"""Independent forcing replay and visible axis-consistency limitation."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plasma_bench.mhd_mesh import MHDMesh
from .electron_reference import ElectronReference
from .current_inertia_reference import g_rate
from .current_thermal_records import provenance
from .verification_records import write_record


def verify(directory=Path('outputs/current_inertia_spatial')):
    data=json.loads((directory/'measurements.json').read_text());assert data['completed']
    assert data['provenance']['engine']==provenance()
    for name,sha in data['provenance']['spatial_sources'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
    assert [c['nr'] for c in data['cases']]==[8,16,32,64,128]
    replay=[];arrays={}
    for case in data['cases']:
        path=directory/case['fields']['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==case['fields']['sha256']
        with np.load(path,allow_pickle=False) as opened:a={k:opened[k] for k in opened.files}
        assert all(np.isfinite(v).all() for v in a.values());arrays[case['nr']]=a
        ref=ElectronReference(**case['reference']);mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],case['nz'])
        np.testing.assert_array_equal(a['radial_faces'],mesh.rf);np.testing.assert_array_equal(a['axial_faces'],mesh.zf)
        # Physical-r quadrature and the full convective E, independent of the
        # measurement's x quadrature and reduced electric field.
        nodes,weights=np.polynomial.legendre.leggauss(8);target=np.zeros_like(a['g'])
        for u,wu in zip(nodes,weights):
            r=mesh.r[:,None]+.5*mesh.dr*u
            for v,wv in zip(nodes,weights):
                z=mesh.z[None]+.5*mesh.dz*v
                target+=.25*wu*wv*r*mesh.dr/mesh.radial_measure[:,None]*g_rate(ref,r,z,case['reduced_mass'],case['delta'],full=True)
        change=float(np.linalg.norm(target-a['reference_g_forcing'])/np.linalg.norm(target))
        assert change<2e-11,(case['nr'],change)
        defect=a['computed_g_forcing']-target
        measured=float(np.sqrt(np.sum(mesh.volumes*defect**2)/np.sum(mesh.volumes*target**2)))
        np.testing.assert_allclose(measured,case['metrics']['g_forcing_relative_L2'],rtol=1e-10)
        replay.append(dict(nr=case['nr'],reference_full_electric_quadrature_relative=change,forcing_relative_L2=measured))
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),constrained_layout=True)
    for key,label in (('g_forcing_relative_L2','forcing, volume L2'),('axis_g_forcing_relative_L2','forcing, first annulus'),('axis_star_correction_relative_L2','curl-w correction, first annulus')):
        axes[0].loglog([c['nr'] for c in data['cases']],[c['metrics'][key] for c in data['cases']],'o-',label=label)
    axes[0].set(title='A volume norm misses the persistent axis error',xlabel='radial cells',ylabel='relative L2 error');axes[0].legend()
    a=arrays[128];z=.5*(a['axial_faces'][:-1]+a['axial_faces'][1:])*1e3
    axes[1].plot(z,a['reference_g_forcing'][0],label='independent continuum')
    axes[1].plot(z,a['computed_g_forcing'][0],'--',label='discrete forcing')
    axes[1].set(title='First-annulus nonlinear forcing, nr=128',xlabel='z [mm]',ylabel='g forcing [T/(m s)]');axes[1].legend()
    fig.suptitle('Non-Beltrami convective inertia: B* forcing before the inertial solve',fontsize=13)
    image=directory/'nonlinear_inertia_axis.png';fig.savefig(image,dpi=160);plt.close(fig)
    names=('examples/verify_current_inertia_spatial.py','tests/test_current_inertia_reference.py')
    report=dict(audit_passed=True,provenance=data['provenance'],replay=replay,
        final_metrics=data['cases'][-1]['metrics'],report_sources={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names},
        figure=dict(file=image.name,sha256=hashlib.sha256(image.read_bytes()).hexdigest()),
        limitation='Persistent first-annulus B* nonlinear-forcing error; actual Bdot after the physical inertia solve must be checked before moving integration.')
    write_record(directory/'report.json',report);print(json.dumps(report['final_metrics'],indent=2),flush=True)


if __name__=='__main__':verify()
