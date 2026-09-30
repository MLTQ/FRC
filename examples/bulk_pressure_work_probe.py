"""Diagnostic conservative pressure stress versus the true kinetic work metric."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.compatible_polynomials import bases
from plasma_bench.compatible_velocity import radial_moments,axial_moments
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from .bulk_thermal_verification import initialize
from .bulk_thermal_records import provenance,sha
from .finite_mass_manufactured import cell_average


def pressure_at(model,species,states,r,z):
    """Single conservative density/continuous entropy trace; no analytic p input."""
    mesh=model.mesh;degree=model.geometry.degree
    rx=radial_moments(mesh,r,degree);rz=axial_moments(mesh,z,degree)
    nx=bases(.5*mesh.rf**2,.5*np.asarray(r)**2,degree)[0]
    nz=bases(mesh.zf,z,degree,periodic=True)[0]
    n=np.asarray(rx@model.density_map.cell_means@rz.T)
    if np.any(n<=0):raise ValueError('Nonpositive diagnostic density trace')
    pressure=np.zeros_like(n)
    for thermal,state in zip(species,states):
        pressure+=n*E*np.exp(np.asarray(nx@state@nz.T))*(n/thermal.density_scale)**(thermal.gamma-1)
    return pressure


def probe(spec):
    ref,model,flux,faces,momentum,species,states=initialize(spec)
    mesh=model.mesh;geo=model.geometry
    if spec.get('uniform_pressure',False):
        r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');n=ref.fields(r,z)['n']
        states=[np.log(10/(n*E*(n/t.density_scale)**(t.gamma-1))) for t in species]
    pr=pressure_at(model,species,states,mesh.rf,geo.axial_points).reshape(mesh.nr+1,mesh.nz,geo.order)
    pz=pressure_at(model,species,states,geo.radial_points,mesh.zf[:-1]).reshape(mesh.nr,geo.order,mesh.nz)
    q,w=np.polynomial.legendre.leggauss(geo.order)
    pr_integrated=2*np.pi*mesh.rf[:,None]*np.sum(pr*(.5*mesh.dz*w),axis=-1)
    pz_integrated=np.sum(pz*geo.wr.reshape(mesh.nr,geo.order,1),axis=1)
    pressure=sum(t.fields(s)['pressure'] for t,s in zip(species,states))
    radial=-(np.diff(pr_integrated,axis=0))/mesh.volumes+model.average(pressure/geo.r)
    axial=-(np.roll(pz_integrated,-1,axis=1)-pz_integrated)/mesh.volumes
    force=np.stack((radial,axial),axis=-1)
    rz=model.density_map.axial;r0=model.density_map.radial;r1=radial_moments(mesh,geo.radial_points,geo.degree,power=1)
    particle=np.stack((geo.tensor(r1,momentum[...,0],rz),geo.tensor(r0,momentum[...,1],rz)),axis=-1)
    velocity=particle/model.density[...,None]
    sampled_force=np.stack((geo.tensor(r1,radial,rz),geo.tensor(r0,axial,rz)),axis=-1)
    mechanical_local=np.sum(velocity*sampled_force,axis=-1)
    rates=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,states)]
    thermal_local=sum(t.fields(s)['internal']*(t.scalar.sample(rate.entropy)+t.gamma*model.density_map.sample(rate.density)/model.density)
                      for t,s,rate in zip(species,states,rates))
    mechanical=float(np.sum(geo.weight*mechanical_local));thermal=float(np.sum(geo.weight*thermal_local))
    scale=float(np.sum(geo.weight*(abs(mechanical_local)+abs(thermal_local))))
    def exact_force(r,z):
        f=ref.fields(r,z);t=ref.thermal(r,z)
        return 2*t['electric'][...,[0,2]]*(E*f['n'])[...,None]
    exact=np.zeros_like(force) if spec.get('uniform_pressure',False) else cell_average(mesh,exact_force,order=18)
    error=np.sum((force-exact)**2,axis=-1);target=np.sum(exact**2,axis=-1)
    relative=float(np.sqrt(np.sum(mesh.volumes*error)/np.sum(mesh.volumes*target))) if np.any(target) else None
    record=dict(spec=spec,conservative_force_relative_L2=relative,
        conservative_force_rms_pa_m=float(np.sqrt(np.sum(mesh.volumes*error)/mesh.volumes.sum())),
        force_axis_rms_pa_m=float(np.sqrt(np.sum(mesh.volumes[0]*error[0])/np.sum(mesh.volumes[0]))),
        mechanical_work=mechanical,thermal_work=thermal,uncancelled_work_scale=scale,work_defect=mechanical+thermal,
        work_backward=abs(mechanical+thermal)/scale,
        net_axial_force=float(np.sum(mesh.volumes*axial)),
        scope='Diagnostic conservative pressure flux only; not adopted as a momentum update or energy correction')
    arrays=dict(density=model.density_map.cell_means,particle_momentum=momentum,electron=states[0],ion=states[1],
        force=force,exact_force=exact,radial_pressure_flux=pr_integrated,axial_pressure_flux=pz_integrated,
        mechanical_cell_work=model.average(mechanical_local),thermal_cell_work=model.average(thermal_local))
    return record,arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);directory=parser.parse_args().directory
    directory.mkdir(parents=True,exist_ok=False);identity=provenance();records=[]
    for index,spec in enumerate([dict(nr=n) for n in (8,16,32,64)]+[dict(nr=32,uniform_pressure=True),dict(nr=32,direction=-1)]):
        record,arrays=probe(spec);path=directory/f'probe_{index:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
        print(record,flush=True)
    source=Path(__file__)
    with tarfile.open(directory/'probe_sources.tar.gz','w:gz') as archive:
        archive.add(source,arcname=source.name);archive.add(source.with_suffix('.md'),arcname=source.with_suffix('.md').name)
    report=dict(calculation_source_sha256=identity['source_sha256'],probe_source_sha256=sha(source),
                probe_companion_sha256=sha(source.with_suffix('.md')),cases=records)
    (directory/'probe.json').write_text(json.dumps(report,indent=2)+'\n')
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Calculation source changed during probe')


if __name__=='__main__':main()
