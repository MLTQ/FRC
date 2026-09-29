"""Measured variable-density induction/energy and graph-pressure consistency."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import traceback
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent,weak_vector
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.compatible_pressure import CompatibleSpeciesPressure
from .finite_mass_manufactured import cell_average
from .variable_current_reference import VariableCurrentReference
from .variable_current_norms import measure
from .variable_current_records import provenance,freeze,sha


def continuity(mesh,reference,order=12):
    """Actual integrated prescribed particle flux, then finite-volume divergence."""
    q,w=np.polynomial.legendre.leggauss(order)
    radial=np.zeros((mesh.nr+1,mesh.nz));axial=np.zeros((mesh.nr,mesh.nz+1))
    for point,weight in zip(q,w):
        r=mesh.rf[:,None];z=mesh.z[None,:]+.5*mesh.dz*point
        radial+=np.pi*mesh.dz*weight*r*reference.fields(r,z)['particle'][...,0]
        r=mesh.r[:,None]+.5*mesh.dr*point;z=mesh.zf[None,:]
        axial+=np.pi*mesh.dr*weight*r*reference.fields(r,z)['particle'][...,2]
    axial[:,-1]=axial[:,0];radial[[0,-1]]=0.
    rate=-(np.diff(radial,axis=0)+np.diff(axial,axis=1))/mesh.volumes
    return rate,radial,axial


def initialize(spec):
    ref=VariableCurrentReference(number_density=spec.get('number_density',1e18),
        contrast=spec.get('contrast',.8),speed=spec.get('speed',200.))
    nr=spec['nr'];mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=12)
    model=VariableCompatibleCurrent(mesh,n,ref.ion_mass,degree=spec.get('degree',5),
                                   quadrature_order=spec.get('quadrature_order',10))
    psi=ref.fields(*mesh.corners)['psi'];psi[[0,-1]]=0.;psi[:,-1]=psi[:,0]
    g=cell_average(mesh,lambda r,z:ref.fields(r,z)['g'],order=12)
    field=HarmonicMagnetic(mesh,psi,g,ref.axial_bias)
    nt,radial,axial=continuity(mesh,ref)
    return ref,model,field,nt,radial,axial


def pressure_data(model,field,reference,current):
    mesh=model.mesh;geo=model.geometry
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');n=reference.fields(r,z)['n']
    arrays={};diagnostics={}
    for label,gamma in (('electron',5/3),('ion',1.4)):
        exact=reference.thermal(r,z,gamma=gamma)
        bridge=CompatibleSpeciesPressure(geo,n,gamma=gamma)
        u=bridge.thermo.from_pressure(exact['pressure']);pair=bridge.pair(field.g,u)
        work=float(current@(model.mass@pair.electric_full))
        thermal=bridge.thermal_energy_rate(pair.entropy_rate)
        old=bridge.thermo.rhs(field.g,u)[1]
        old_work=bridge.thermal_energy_rate(old)
        volume=bridge.geometry.volumes
        scale=float(abs(current)@abs(model.mass@pair.electric_full))
        entropy=volume*n/u*pair.entropy_rate
        diagnostics[label]=dict(work_w=work,thermal_rate_w=thermal,
            work_backward=abs(work-thermal)/max(scale,abs(thermal),np.finfo(float).tiny),
            wrong_old_stream_work_relative=abs(old_work-work)/max(abs(thermal),np.finfo(float).tiny),
            entropy_backward=float(abs(entropy.sum())/max(np.sum(abs(entropy)),np.finfo(float).tiny)),
            nodal_entropy_rate_relative_L2=float(np.sqrt(np.sum(volume*(pair.entropy_rate-exact['entropy_rate'])**2)/np.sum(volume*exact['entropy_rate']**2))))
        arrays[label+'_electric']=pair.electric_full
        arrays[label+'_entropy_rate']=pair.entropy_rate
    return arrays,diagnostics


def run_case(spec):
    ref,m,field,nt,radial,axial=initialize(spec);geo=m.geometry
    current=field.current(m)
    other=weak_vector(geo,ref.fields(geo.r,geo.z)['other_forcing'])
    density=m.density_forcing(current,nt)
    solved=KernelInduction(m,tolerance=1e-10,max_iterations=3000).solve(other+density)
    arrays=dict(magnetic=field.packed_increment(geo),current=current,magnetic_rate=solved.magnetic,
        current_rate=solved.current,electric=solved.electric,n_cell=m.density_map.cell_means,
        nt_cell=nt,particle_radial=radial,particle_axial=axial,forcing_other=other,forcing_density=density)
    pressure,pressure_diagnostics=pressure_data(m,field,ref,current);arrays.update(pressure)
    krate=m.current_energy_rate(current,solved.current,nt)
    production=dict(current_energy=m.current_energy(current),current_energy_rate=krate['total'],
        density_energy_rate=krate['density'],magnetic_energy_rate=m.magnetic_energy_rate(field,solved.magnetic))
    arrays.update({'production_'+k:v for k,v in production.items()})
    norms,independent=measure(m.mesh,ref,geo.degree,arrays,order=12)
    higher,higher_cells=measure(m.mesh,ref,geo.degree,arrays,order=16)
    cell_audit={k:float(np.linalg.norm(production[k]-independent[k])/max(np.linalg.norm(independent[k]),np.finfo(float).tiny)) for k in production}
    quadrature={k:float(np.linalg.norm(higher_cells[k]-independent[k])/max(np.linalg.norm(higher_cells[k]),np.finfo(float).tiny)) for k in production}
    volume=m.mesh.volumes
    total_rate=float(np.sum(volume*(production['current_energy_rate']+production['magnetic_energy_rate'])))
    target=float(-current@other-np.sum(volume*krate['density']))
    energy_scale=float(abs(current)@abs(other)+np.sum(volume*(abs(production['magnetic_energy_rate'])+abs(production['current_energy_rate']))))
    divergence=-(np.diff(radial,axis=0)+np.diff(axial,axis=1))
    nt_exact=cell_average(m.mesh,lambda r,z:ref.fields(r,z)['nt'],order=12)
    diagnostics=dict(spec=spec,reference=asdict(ref),status='complete',norms=norms,
        higher_quadrature_norms=higher,pressure=pressure_diagnostics,
        production_energy_relative_difference=cell_audit,independent_quadrature_relative_difference=quadrature,
        energy_exchange_backward=abs(total_rate-target)/max(energy_scale,np.finfo(float).tiny),
        energy_rate_w=total_rate,energy_exchange_target_w=target,
        continuity_relative_L2=float(np.linalg.norm(nt-nt_exact)/np.linalg.norm(nt_exact)),
        continuity_integral_backward=float(abs(divergence.sum())/np.sum(abs(divergence))),
        residual=solved.residual,residual_components=solved.residual_components,
        iterations=solved.iterations,refinements=solved.refinements,
        min_density_m3=float(m.density.min()),max_density_m3=float(m.density.max()))
    return diagnostics,arrays


def protocol():
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32,64)]
    specs.extend(dict(nr=32,quadrature_order=q) for q in (8,12))
    specs.extend(dict(nr=32,number_density=n) for n in (1e16,1e20))
    specs.extend([dict(nr=32,contrast=0.),dict(nr=32,speed=-200.)])
    return specs


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    args=parser.parse_args();directory=args.directory;directory.mkdir(parents=True,exist_ok=False)
    identity=provenance();freeze(directory,identity);records=[]
    for index,spec in enumerate(protocol()):
        name=f'case_{index:02d}';print(name,spec,flush=True)
        try:
            result,arrays=run_case(spec)
            path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
            result.update(arrays_file=path.name,arrays_sha256=sha(path))
        except Exception as error:
            result=dict(spec=spec,status='failed',error=repr(error),traceback=traceback.format_exc())
        result['source_sha256']=identity['source_sha256'];records.append(result)
        (directory/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        print(result['status'],result.get('norms',{}).get('metrics',{}).get('magnetic_rate',{}),flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')
    (directory/'records.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['status']!='complete' for r in records):raise RuntimeError('Preserved failed cases; study is incomplete')


if __name__=='__main__':main()
