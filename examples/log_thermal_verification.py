"""Full physical thermal rate protocol with density work and original induction."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import traceback
import numpy as np
from plasma_bench.compatible_pressure import CompatibleSpeciesPressure
from plasma_bench.kernel_induction import KernelInduction
from .log_thermal_cases import initialize,protocol,target
from .log_thermal_norms import measure
from .log_thermal_records import provenance,freeze,sha


def relative(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),np.finfo(float).tiny))


def run_case(spec):
    ref,m,field,nt,radial,axial,species,states=initialize(spec);geo=m.geometry;mesh=m.mesh
    current=field.current(m);pairs={label:t.pair(current,states[label]) for label,t in species.items()}
    arrays=dict(magnetic=field.packed_increment(geo),current=current,n_cell=m.density_map.cell_means,
                nt_cell=nt,particle_radial=radial,particle_axial=axial)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');nodal_n=ref.fields(r,z)['n']
    production={};work={}
    for label,t in species.items():
        s=states[label];pair=pairs[label];kind=spec.get('kind','smooth')
        legacy=CompatibleSpeciesPressure(geo,nodal_n,gamma=t.gamma)
        old=legacy.thermo.from_pressure(target(ref,r,z,t.gamma,kind)['pressure'])
        old_pair=legacy.pair(field.g,old)
        arrays.update({label+'_entropy':s,label+'_rate':pair.entropy_rate,
            label+'_electric':pair.electric_full,label+'_reduced':pair.electric_reduced,
            label+'_potential':pair.potential,label+'_energy_per_particle':pair.energy_per_particle,
            label+'_legacy_rate':old_pair.entropy_rate,label+'_legacy_proxy':old})
        rates=t.energy_rate(s,pair.entropy_rate,nt)
        production[label]=dict(internal=t.energy_density(s),internal_entropy_rate=rates['entropy'],
                               internal_density_rate=rates['density'])
        power=float(current@(m.mass@pair.electric_full));thermal=float(np.sum(mesh.volumes*rates['entropy']))
        scale=float(abs(current)@abs(m.mass@pair.electric_full)+np.sum(mesh.volumes*abs(rates['entropy'])))
        entropy=t.entropy_rate(s,pair.entropy_rate,np.zeros_like(nt))
        proxy=t.proxy_rate(s,pair.entropy_rate,np.zeros_like(nt))
        work[label]=dict(pressure_power_w=power,internal_entropy_rate_w=thermal,
            pressure_work_backward=abs(power-thermal)/max(scale,np.finfo(float).tiny),
            physical_entropy_backward=abs(float(np.sum(mesh.volumes*entropy)))/max(float(np.sum(mesh.volumes*abs(entropy))),np.finfo(float).tiny),
            proxy_integral_defect=float(np.sum(mesh.volumes*proxy)),
            proxy_defect_relative=float(np.sum(mesh.volumes*proxy))/max(float(np.sum(mesh.volumes*abs(proxy))),np.finfo(float).tiny),
            scalar_residuals=pair.residuals)
    reduced=m.beta*pairs['electron'].electric_reduced-m.alpha*pairs['ion'].electric_reduced
    potential=m.beta*pairs['electron'].potential-m.alpha*pairs['ion'].potential
    forcing=m.mass@reduced+m.density_forcing(current,nt)
    solved=KernelInduction(m,tolerance=1e-10,max_iterations=3000).solve(forcing)
    full=solved.electric+geo.gradient(potential)
    arrays.update(joint_magnetic_rate=solved.magnetic,joint_current_rate=solved.current,
        joint_electric=full,joint_forcing_full=forcing+m.mass@geo.gradient(potential),
        joint_electric_reduced=solved.electric,joint_forcing_reduced=forcing)
    kinetic=m.current_energy_rate(current,solved.current,nt)
    production['field']=dict(current_energy_rate=kinetic['total'],density_energy_rate=kinetic['density'],
                             magnetic_energy_rate=m.magnetic_energy_rate(field,solved.magnetic))
    for label,values in production.items():
        arrays.update({'production_'+label+'_'+key:value for key,value in values.items()})
    thermal=m.beta*production['electron']['internal_entropy_rate']-m.alpha*production['ion']['internal_entropy_rate']
    density=production['electron']['internal_density_rate']+production['ion']['internal_density_rate']
    total=production['field']['magnetic_energy_rate']+kinetic['total']+thermal+density
    wanted=-kinetic['density']+density
    scale=np.sum(mesh.volumes*(abs(production['field']['magnetic_energy_rate'])+abs(kinetic['total'])+abs(thermal)+abs(density)))
    exchange=float(abs(np.sum(mesh.volumes*(total-wanted)))/max(scale,np.finfo(float).tiny))
    norms,independent=measure(mesh,ref,spec,arrays,order=12)
    higher,higher_cells=measure(mesh,ref,spec,arrays,order=16)
    energy_audit={label:{key:relative(production[label][key]-value,value) for key,value in values.items()} for label,values in independent.items()}
    quadrature={label:{key:relative(higher_cells[label][key]-value,higher_cells[label][key]) for key,value in values.items()} for label,values in independent.items()}
    record=dict(status='complete',spec=spec,reference=asdict(ref),norms=norms,higher_quadrature_norms=higher,
        work=work,energy_exchange_backward=exchange,independent_energy_agreement=energy_audit,
        independent_quadrature_agreement=quadrature,induction_residual=solved.residual,
        induction_components=solved.residual_components,induction_iterations=solved.iterations)
    return record,arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory;directory.mkdir(parents=True,exist_ok=False)
    identity=provenance();freeze(directory,identity);records=[]
    for index,spec in enumerate(protocol()):
        name=f'rate_{index:02d}';print(name,spec,flush=True)
        try:
            record,arrays=run_case(spec);path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
            record.update(arrays_file=path.name,arrays_sha256=sha(path))
        except Exception as error:
            record=dict(status='failed',spec=spec,error=repr(error),traceback=traceback.format_exc())
        record['source_sha256']=identity['source_sha256'];records.append(record)
        (directory/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n')
        print(record['status'],record.get('norms',{}).get('nodal',{}).get('electron',{}).get('physical',{}).get('axis_four',{}),flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during verification')
    (directory/'rates.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['status']!='complete' for r in records):raise RuntimeError('Preserved failed rate cases')


if __name__=='__main__':main()
