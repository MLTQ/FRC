"""Preserve the distinction between operator cancellation and rounded input loss."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.compatible_current import CompatibleCurrent
from plasma_bench.kernel_induction import KernelInduction
from .compatible_wave_verification import initial
from .compatible_wave_records import provenance
from .verification_records import write_record


def run(directory='outputs/compatible_wave_verification'):
    directory=Path(directory);previous=Path('outputs/compatible_wave_rounded_attempt')
    row=json.loads((previous/'branch0_32.json').read_text());ref,mesh,_,_,g=initial(row['configuration'])
    model=CompatibleCurrent(mesh,ref.number_density,ref.ion_mass,degree=5,quadrature_order=10);geo=model.geometry
    with np.load(previous/row['fields']['file']) as data:total=data['magnetic'];v=data['velocity']
    assert hashlib.sha256((previous/row['fields']['file']).read_bytes()).hexdigest()==row['fields']['sha256']
    harmonic=geo.pack(.5*ref.bias*mesh.corners[0]**2,np.zeros_like(g))
    psi=ref.fields(*mesh.corners)['psi'];psi[0]=0.;psi[:,-1]=psi[:,0]
    target=ref.derivatives(geo.r,geo.z)[2]['magnetic'].real
    norm=lambda value:float(np.sqrt(np.sum(geo.weight[...,None]*value*value)))
    rows={};arrays=dict(total=total,harmonic=harmonic,velocity=v)
    for name,field in (('rounded_total',total),('subtract_after_rounding',total-harmonic),('separate_increment',geo.pack(psi,g))):
        current=model.current(field)
        force=model.mass@geo.project_electric(np.cross(-v+model.delta*geo.sample_electric(current)/(ref.number_density*E),[0.,0.,ref.bias]))
        solved=KernelInduction(model,tolerance=1e-10,max_iterations=2000).solve(force)
        actual=geo.magnetic_samples(*geo.unpack(solved.magnetic))
        rows[name]=dict(magnetic_rate_relative_L2=norm(actual-target)/norm(target),linear_residual=solved.residual)
        arrays[name+'_input']=field;arrays[name+'_current']=current;arrays[name+'_rate']=solved.magnetic
    j0=model.current(harmonic);arrays['roundoff_harmonic_current']=j0
    path=directory/'bias_diagnosis.npz';np.savez_compressed(path,**arrays)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    record=dict(scope='Diagnostic only: same rounded input versus separately retained increment. No physical state or velocity is fitted.',
        reference_configuration=row['configuration'],input_record_sha256=sha(previous/'branch0_32.json'),input_field_sha256=row['fields']['sha256'],
        current_source=provenance(),helper_sources={p:sha(p) for p in ('examples/compatible_bias_diagnosis.py','examples/compatible_bias_diagnosis.md')},
        cases=rows,roundoff_harmonic_current_rms_a_m2=norm(geo.sample_electric(j0))/np.sqrt(mesh.volumes.sum()),
        fields=dict(file=path.name,sha256=sha(path)))
    write_record(directory/'bias_diagnosis.json',record);print(json.dumps(rows,indent=2))


if __name__=='__main__':run()
