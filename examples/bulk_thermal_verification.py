"""Matched conservative-continuity and physical bulk thermal rate protocol."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import traceback
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from .variable_current_reference import VariableCurrentReference
from .finite_mass_manufactured import cell_average
from .bulk_thermal_norms import rate_norms
from .bulk_thermal_records import provenance,freeze,sha


def initialize(spec):
    ref=VariableCurrentReference(number_density=spec.get('number_density',1e18),speed=200*spec.get('direction',1),contrast=spec.get('contrast',.8))
    nr=spec['nr'];mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=14)
    model=VariableCompatibleCurrent(mesh,n,ref.ion_mass,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    flux=CompatibleParticleFlux(model.density_map)
    momentum=cell_average(mesh,lambda r,z:ref.fields(r,z)['particle'][...,[0,2]],order=14)
    faces=flux.from_momentum(momentum);r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij')
    nnode=ref.fields(r,z)['n'];species=[];states=[]
    for gamma in (5/3,1.4):
        thermal=CompatibleLogThermal(model,gamma=gamma,scalar=species[0].scalar if species else None)
        pressure=ref.thermal(r,z,gamma=gamma)['pressure']
        s=np.full_like(nnode,2.) if spec.get('constant_entropy',False) else np.log(pressure/(nnode*E*(nnode/1e20)**(gamma-1)))
        species.append(thermal);states.append(s)
    return ref,model,flux,faces,momentum,species,states


def calculate(spec):
    ref,model,flux,faces,momentum,species,states=initialize(spec)
    arrays=dict(density=model.density_map.cell_means,particle_momentum=momentum,radial=faces.radial,axial=faces.axial,density_rate=flux.density_rate(faces))
    work={}
    for label,thermal,state in zip(('electron','ion'),species,states):
        rate=CompatibleBulkThermal(thermal,flux).rate(state,faces)
        energy=thermal.energy_rate(state,rate.entropy,rate.density)
        thermal_rate=float(np.sum(model.mesh.volumes*energy['total']))
        face_rate=float(np.sum(faces.radial*rate.face_work.radial)+np.sum(faces.axial*rate.face_work.axial))
        scale=float(np.sum(model.mesh.volumes*(abs(energy['entropy'])+abs(energy['density']))))
        work[label]=dict(thermal_rate=thermal_rate,face_rate=face_rate,backward=abs(thermal_rate-face_rate)/scale,residuals=rate.residuals)
        arrays.update({label:state,label+'_rate':rate.entropy,label+'_energy_rate':energy['total'],
                       label+'_radial_work':rate.face_work.radial,label+'_axial_work':rate.face_work.axial})
    norms,local=rate_norms(model.mesh,ref,spec,arrays,order=14)
    higher,_=rate_norms(model.mesh,ref,spec,arrays,order=16)
    agreement={label:float(np.linalg.norm(local[label]-arrays[label+'_energy_rate'])/np.linalg.norm(local[label])) for label in ('electron','ion')}
    return dict(status='complete',spec=spec,reference=asdict(ref),work=work,norms=norms,higher_quadrature_norms=higher,independent_local_energy=agreement),arrays


def protocol():
    return ([dict(nr=n,degree=d) for d in (3,5) for n in (8,16,32,64)]
        +[dict(nr=32,quadrature_order=16),dict(nr=32,number_density=1e16),dict(nr=32,number_density=1e20),
          dict(nr=32,direction=-1),dict(nr=32,constant_entropy=True),dict(nr=32,contrast=1.4)])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);directory=parser.parse_args().directory
    directory.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(directory,identity);records=[]
    for index,spec in enumerate(protocol()):
        name=f'rate_{index:02d}';print(name,spec,flush=True)
        try:
            record,arrays=calculate(spec);path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
            record.update(arrays_file=path.name,arrays_sha256=sha(path))
        except Exception as error:record=dict(status='failed',spec=spec,error=repr(error),traceback=traceback.format_exc())
        record['source_sha256']=identity['source_sha256'];records.append(record)
        (directory/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(record['status'],flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during rates')
    (directory/'rates.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['status']!='complete' for r in records):raise RuntimeError('Preserved failed rate cases')


if __name__=='__main__':main()
