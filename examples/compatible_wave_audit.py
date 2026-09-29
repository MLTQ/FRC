"""Read-only original-equation and physical-norm audit of compatible wave records."""
from pathlib import Path
import hashlib
import json
import tarfile
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.compatible_current import CompatibleCurrent
from plasma_bench.compatible_velocity import velocity
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .compatible_wave_verification import initial
from .compatible_field_norms import measure


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(directory,*,current):
    directory=Path(directory);record=json.loads((directory/'source.json').read_text())
    assert digest(directory/'source_snapshot.tar.gz')==record['archive_sha256']
    changed=[]
    with tarfile.open(directory/'source_snapshot.tar.gz') as tar:
        for name,sha in record['provenance']['source_files'].items():
            assert hashlib.sha256(tar.extractfile(name).read()).hexdigest()==sha,name
            if digest(name)!=sha:changed.append(name)
    if current:assert not changed,changed
    else:assert set(changed)<= {'examples/compatible_wave_verification.py','examples/compatible_field_norms.py'},changed
    return record['provenance'],changed


def case(directory,name,frozen):
    directory=Path(directory);record=json.loads((directory/(name+'.json')).read_text())
    assert record['provenance']==frozen and record['completed']
    assert digest(directory/record['fields']['file'])==record['fields']['sha256']
    with np.load(directory/record['fields']['file']) as a:arrays={key:a[key] for key in a.files}
    configuration=record['configuration'];ref,mesh,matter,psi,g=initial(configuration)
    model=CompatibleCurrent(mesh,ref.number_density,ref.ion_mass,degree=configuration['degree'],quadrature_order=configuration['quadrature_order'])
    geo=model.geometry;magnetic=geo.pack(psi,g)
    np.testing.assert_array_equal(matter,arrays['matter']);np.testing.assert_array_equal(magnetic,arrays['magnetic'])
    np.testing.assert_array_equal(model.current(magnetic),arrays['current'])
    v=velocity(geo,matter,degree=configuration['degree'])
    np.testing.assert_array_equal(v,arrays['velocity'])
    force=model.mass@geo.project_electric(np.cross(-v+model.delta*geo.sample_electric(arrays['current'])/(ref.number_density*E),[0.,0.,ref.bias]))
    np.testing.assert_array_equal(force,arrays['forcing'])
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');reference=ref.fields(r,z)
    potential=(model.alpha*reference['ion_pressure']-model.beta*reference['electron_pressure'])/(ref.number_density*E)
    np.testing.assert_array_equal(potential,arrays['potential'])
    energy=.5*magnetic@(model.hessian@magnetic)
    if configuration.get('representation')=='affine_harmonic':
        field=HarmonicMagnetic(mesh,psi,g,ref.bias if configuration['bias_in_state'] else 0.)
        assert field.axial_bias==arrays['axial_bias_t']
        np.testing.assert_array_equal(field.rounded_total(geo),arrays['rounded_total'])
        assert field.energy_components(model)==record['magnetic_energy_components_j']
        energy=sum(field.energy_components(model).values())
    result=measure(mesh,ref,configuration['degree'],arrays,bias_in_state=configuration['bias_in_state'],order=14)
    for key in result['metrics']:
        for norm,value in result['metrics'][key].items():np.testing.assert_allclose(value,record['metrics'][key][norm],rtol=1e-13,atol=1e-30)
    magnetic_energy_error=abs(result['magnetic_energy_j']-energy)/abs(energy)
    current_energy=.5*arrays['current']@(model.inertia@arrays['current'])
    current_energy_error=abs(result['current_energy_j']-current_energy)/abs(current_energy)
    increment_energy_error=None
    if configuration.get('representation')=='affine_harmonic':
        increment=measure(mesh,ref,configuration['degree'],dict(arrays,axial_bias_t=np.array(0.)),bias_in_state=configuration['bias_in_state'],order=14)
        expected=field.energy_components(model)['increment']
        increment_energy_error=abs(increment['magnetic_energy_j']-expected)/abs(expected)
        assert increment_energy_error<1e-8
    # Physical norms at a different order need not be machine-identical to
    # assembled mass quadrature. Their measured differences are explicit.
    assert magnetic_energy_error<1e-8 and current_energy_error<1e-8
    rate=arrays['magnetic_rate'];jrate=arrays['current_rate'];electric=arrays['electric_full']
    norm=lambda x:float(np.linalg.norm(x))
    gradient=geo.gradient(arrays['potential'])
    ohm=model.mass@electric-model.inertia@jrate-arrays['forcing']-model.mass@gradient
    ohm_scale=abs(model.mass@electric)+abs(model.inertia@jrate)+abs(arrays['forcing'])+abs(model.mass@gradient)
    ampere=model.mass@jrate+geo.faraday.T@(model.hessian@rate)
    ampere_scale=abs(model.mass@jrate)+abs(geo.faraday.T@(model.hessian@rate))
    faraday=geo.faraday@electric-rate
    faraday_scale=abs(geo.faraday)@abs(electric)+abs(rate)
    backward=dict(ohm=norm(ohm)/norm(ohm_scale),ampere_rate=norm(ampere)/norm(ampere_scale),faraday=norm(faraday)/norm(faraday_scale))
    assert max(backward.values())<1e-10
    full_rate_relative=float(np.sqrt(abs(faraday@(model.hessian@faraday))/abs(rate@(model.hessian@rate))))
    assert record['linear_residual']<=1e-10
    np.testing.assert_allclose(full_rate_relative,record['original_full_faraday_rate_relative'],rtol=1e-13)
    return record,dict(record_sha256=digest(directory/(name+'.json')),backward=backward,
        full_faraday_rate_relative=full_rate_relative,magnetic_energy_relative_difference=magnetic_energy_error,
        current_energy_relative_difference=current_energy_error,increment_energy_relative_difference=increment_energy_error)
