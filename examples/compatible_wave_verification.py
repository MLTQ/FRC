"""Physical wave-rate convergence for high-order compatible magnetic/current fields."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import tarfile
import numpy as np
from scipy.constants import elementary_charge as E,atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_current import CompatibleCurrent
from plasma_bench.compatible_velocity import velocity
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_manufactured import cell_average
from .compatible_field_norms import measure
from .compatible_wave_records import provenance
from .verification_records import write_record


def protocol():
    def case(name,nr,branch=0,degree=5,order=10,amplitude=1.,bias=True,density=1e18,sense=1):
        ref=FiniteMassWave(.035,.24,density,39.948*atomic_mass,.005,6*density*E,gamma_i=1.4,branch=branch,amplitude=amplitude,sense=sense)
        return dict(name=name,nr=nr,nz=2*nr,degree=degree,quadrature_order=order,bias_in_state=bias,reference=asdict(ref))
    cases=[case(f'branch{branch}_{nr}',nr,branch) for branch in range(3) for nr in (8,16,32,64)]
    cases += [case(f'cubic_slow_{nr}',nr,degree=3) for nr in (8,16,32,64)]
    cases += [case(f'matched_linear_{nr}',nr,amplitude=.01,bias=False) for nr in (8,16,32,64)]
    cases += [case('quadrature8',32,order=8),case('quadrature12',32,order=12),case('weak_full_bias',32,amplitude=.01),case('original_density',32,density=1e20),case('opposite_slow',32,sense=-1)]
    return [dict(case,representation='affine_harmonic') for case in cases]


def initial(case):
    ref=FiniteMassWave(**case['reference']);mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],case['nz'])
    def momenta(r,z):
        v=ref.fields(r,z)['velocity'];rho=np.full_like(r,ref.density)
        return np.stack((rho,rho*v[...,0],rho*v[...,2],np.zeros_like(r),rho*r*v[...,1]),axis=-1)
    matter=cell_average(mesh,momenta,order=12)
    r,z=mesh.corners;psi=ref.fields(r,z)['psi']
    if case.get('representation')=='affine_harmonic':psi[-1]=0.
    elif case['bias_in_state']:psi+=.5*ref.bias*r*r
    psi[0]=0.;psi[:,-1]=psi[:,0]
    g=cell_average(mesh,lambda r,z:ref.fields(r,z)['magnetic'][...,1]/r,order=12)
    return ref,mesh,matter,psi,g


def run_case(case):
    ref,mesh,matter,psi,g=initial(case)
    model=CompatibleCurrent(mesh,ref.number_density,ref.ion_mass,degree=case['degree'],quadrature_order=case['quadrature_order'])
    geo=model.geometry;magnetic=geo.pack(psi,g)
    field=HarmonicMagnetic(mesh,psi,g,ref.bias if case['bias_in_state'] else 0.) if case.get('representation')=='affine_harmonic' else None
    current=field.current(model) if field is not None else model.current(magnetic)
    v=velocity(geo,matter,degree=case['degree']);j=geo.sample_electric(current)
    force=model.mass@geo.project_electric(np.cross(-v+model.delta*j/(ref.number_density*E),[0.,0.,ref.bias]))
    result=KernelInduction(model,tolerance=1e-10,max_iterations=2000).solve(force)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');fields=ref.fields(r,z)
    potential=(model.alpha*fields['ion_pressure']-model.beta*fields['electron_pressure'])/(ref.number_density*E)
    gradient=geo.gradient(potential);electric=result.electric+gradient
    arrays=dict(matter=matter,magnetic=magnetic,current=current,magnetic_rate=result.magnetic,current_rate=result.current,
        electric_reduced=result.electric,electric_full=electric,forcing=force,potential=potential,velocity=v)
    if field is not None:arrays.update(axial_bias_t=np.array(field.axial_bias),rounded_total=field.rounded_total(geo))
    a=measure(mesh,ref,case['degree'],arrays,bias_in_state=case['bias_in_state'],order=10)
    b=measure(mesh,ref,case['degree'],arrays,bias_in_state=case['bias_in_state'],order=14)
    q=max(abs(a['metrics'][key]['relative_L2']-b['metrics'][key]['relative_L2']) for key in a['metrics'])
    norm=lambda value:float(np.linalg.norm(value))
    ohm=model.mass@electric-model.inertia@result.current-force-model.mass@gradient
    backward=norm(ohm)/norm(abs(model.mass@electric)+abs(model.inertia@result.current)+abs(force)+abs(model.mass@gradient))
    curl=geo.faraday@electric-result.magnetic
    curl_backward=norm(curl)/norm(abs(geo.faraday)@abs(electric)+abs(result.magnetic))
    curl_physical=float(np.sqrt(abs(curl@(model.hessian@curl))/abs(result.magnetic@(model.hessian@result.magnetic))))
    physical_norm=lambda value:float(np.sqrt(np.sum(geo.weight[...,None]*value*value)))
    exact=ref.fields(geo.r,geo.z)
    energy_components=field.energy_components(model) if field is not None else dict(total=float(.5*magnetic@(model.hessian@magnetic)))
    magnetic_energy=sum(energy_components.values())
    current_energy=float(.5*current@(model.inertia@current))
    total=arrays.get('rounded_total',magnetic)
    ampere=model.mass@current+geo.faraday.T@(model.hessian@total)
    ampere_scale=abs(model.mass)@abs(current)+abs(geo.faraday.T)@(abs(model.hessian)@abs(total))
    record=dict(completed=True,configuration=case,metrics=b['metrics'],quadrature_absolute_change=q,
        velocity_relative_L2=physical_norm(v-exact['velocity'])/physical_norm(exact['velocity']),
        original_full_ohm_backward=backward,original_full_faraday_backward=curl_backward,original_full_faraday_rate_relative=curl_physical,
        linear_residual=result.residual,linear_iterations=result.iterations,linear_refinements=result.refinements,induction_components=result.residual_components,
        physical_magnetic_energy_j=b['magnetic_energy_j'],matrix_magnetic_energy_j=magnetic_energy,
        physical_current_energy_j=b['current_energy_j'],matrix_current_energy_j=current_energy,magnetic_energy_components_j=energy_components,
        original_full_ampere_backward=norm(ampere)/norm(ampere_scale),rounded_ampere_current_relative=norm(ampere)/norm(model.mass@current))
    assert backward<1e-10 and curl_backward<1e-10 and record['original_full_ampere_backward']<1e-10
    return record,arrays


def run(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);frozen=provenance()
    assert not (directory/'source.json').exists(),'Use a fresh directory; never overwrite an earlier attempt'
    snapshot=directory/'source_snapshot.tar.gz'
    with tarfile.open(snapshot,'w:gz') as tar:
        for name in frozen['source_files']:tar.add(name,arcname=name)
    write_record(directory/'source.json',dict(provenance=frozen,archive_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest()))
    for case in protocol():
        start=perf_counter()
        # Any exception is terminal and visible; completed records are retained.
        # Do not relabel a reporting failure as a physical initialization failure.
        record,arrays=run_case(case);assert provenance()==frozen
        fields=directory/(case['name']+'.npz');np.savez_compressed(fields,**arrays)
        record.update(provenance=frozen,elapsed_s=perf_counter()-start,fields=dict(file=fields.name,sha256=hashlib.sha256(fields.read_bytes()).hexdigest()))
        write_record(directory/(case['name']+'.json'),record)
        print(case['name'],record['metrics']['magnetic_rate'],record['quadrature_absolute_change'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/compatible_wave_verification')
    run(parser.parse_args().directory)
