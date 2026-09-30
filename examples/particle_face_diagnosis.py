"""Independent all-mode and full-transit diagnosis of normal-face bias."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_velocity import axial_moments
from plasma_bench.normal_momentum_faces import axial_faces
from plasma_bench.particle_flux_basis import difference
from .bulk_thermal_records import provenance,sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);directory=parser.parse_args().directory
    directory.mkdir(parents=True,exist_ok=False);records=[]
    for degree in (3,5):
        for nz in (16,32,64):
            mesh=MHDMesh(.035,.12,nz//2,nz);k=np.pi/mesh.half_length
            initial=1+.25*np.diff(np.sin(k*mesh.zf))/(k*mesh.dz)
            for direction in (1,-1):
                speed=2000*direction;period=2*mesh.half_length/abs(speed);steps=5*nz;dt=period/steps
                for kind,face in (('initial',axial_moments(mesh,mesh.zf[:-1],degree)),('centered',axial_faces(mesh,degree))):
                    operator=(-speed*difference(mesh.zf,True)@face).toarray();values=initial.copy();snapshots=[values.copy()]
                    for step in range(1,steps+1):
                        a=operator@values;b=operator@(values+.5*dt*a);c=operator@(values+.5*dt*b);d=operator@(values+dt*c)
                        values=values+dt*(a+2*b+2*c+d)/6
                        if step%(steps//4)==0:snapshots.append(values.copy())
                    eigen=np.linalg.eigvals(operator)
                    record=dict(kind=kind,degree=degree,nz=nz,direction=direction,
                        growth_per_transit=float(eigen.real.max()*period),
                        density_wave_error=float(np.linalg.norm(values-initial)/np.linalg.norm(initial-1)),
                        min_density=float(values.min()),max_density=float(values.max()),
                        particle_drift=float(abs(values.sum()/initial.sum()-1)),
                        symmetric_part_relative=float(np.linalg.norm(operator+operator.T)/np.linalg.norm(operator)))
                    path=directory/f'case_{len(records):02d}.npz'
                    np.savez_compressed(path,operator=operator,initial=initial,snapshots=np.stack(snapshots),eigenvalues=eigen)
                    record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
    root=Path(__file__).resolve().parents[1];helper=Path(__file__)
    paths=[helper,helper.with_suffix('.md'),root/'plasma_bench/normal_momentum_faces.py',root/'plasma_bench/normal_momentum_faces.md']
    with tarfile.open(directory/'diagnosis_sources.tar.gz','w:gz') as archive:
        for path in paths:archive.add(path,arcname=str(path.relative_to(root)))
    result=dict(bulk_calculation_source_sha256=provenance()['source_sha256'],additional_source_sha256={str(p.relative_to(root)):sha(p) for p in paths},cases=records,
        scope='Independent linear particle-advection diagnosis; neither physical positivity nor full nonlinear stability is asserted')
    (directory/'diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in records:print({k:v for k,v in r.items() if k not in ('arrays_file','arrays_sha256')})


if __name__=='__main__':main()
