"""Independent polynomial test of radial velocity/continuity metric consistency."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.compatible_density import CompatibleDensity
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux
from plasma_bench.conservative_kinematics import ConservativeKinematics
from .radial_work_records import provenance,sha


def calculate(nr,degree):
    mesh=MHDMesh(.035,.12,nr,2*nr);g=CompatibleGeometry(mesh,degree=degree,quadrature_order=16)
    density=CompatibleDensity(g,np.ones((nr,2*nr)));flux=CompatibleParticleFlux(density)
    kin=ConservativeKinematics(g,1.);r=g.radial_points;w=g.wr
    R=kin.radial[1].toarray();Q=flux.hx[:,1:-1].toarray()/r[:,None]
    B=(mesh.rf[:,None]*flux.radial_face.toarray())[1:-1];L=Q@B
    HR=R.T@(w[:,None]*R);HQ=Q.T@(w[:,None]*Q);HL=L.T@(w[:,None]*L)
    # Physical force is regular at the axis and vanishes at the impermeable wall.
    # Its primitive pressure is an even polynomial. No PDE solver supplies it.
    target=r/mesh.radius*(1-(r/mesh.radius)**2)**2
    cell=np.sum((w*target).reshape(nr,g.order),axis=1)/(2*np.pi*mesh.radial_measure)
    actual_cell=solve(HR,L.T@(w*target),assume_a='pos')
    actual_face=solve(HQ,Q.T@(w*target),assume_a='pos')
    sampled_cell=R@actual_cell;sampled_face=Q@actual_face
    def norm(value):
        return dict(relative_L2=float(np.sqrt(np.sum(w*(value-target)**2)/np.sum(w*target**2))),
                    axis_relative_L2=float(np.sqrt(np.sum(w[:g.order]*(value[:g.order]-target[:g.order])**2)/np.sum(w[:g.order]*target[:g.order]**2))),
                    maximum=float(np.max(abs(value-target))))
    records=dict(nr=nr,degree=degree,cell_metric_force=norm(sampled_cell),face_metric_force=norm(sampled_face),
                 direct_cell_representation=norm(R@cell),continuity_flux_representation=norm(L@cell),
                 moment_to_face_rank=int(np.linalg.matrix_rank(B)),cell_flux_metric_eigenvalues=np.linalg.eigvalsh(HL).tolist(),
                 face_metric_condition=float(np.linalg.cond(HQ)))
    arrays=dict(radius=r,weight=w,target=target,cell_momentum_basis=R,face_momentum_basis=Q,cell_to_face=B,
                cell_metric=HR,face_metric=HQ,flux_cell_metric=HL,cell_force=actual_cell,face_force=actual_face,
                sampled_cell=sampled_cell,sampled_face=sampled_face,analytic_cell_means=cell)
    return records,arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);records=[];source=Path(__file__);identity=provenance()
    for degree in (3,5):
        for nr in (8,16,32,64):
            record,arrays=calculate(nr,degree);path=d/f'metric_d{degree}_{nr}.npz';np.savez_compressed(path,**arrays)
            record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
            print(degree,nr,record['cell_metric_force'],record['face_metric_force'],flush=True)
    (d/'records.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],probe_sha256=sha(source),records=records),indent=2)+'\n')
    (d/'probe_source.py').write_bytes(source.read_bytes());(d/'probe_source.md').write_bytes(source.with_suffix('.md').read_bytes())


if __name__=='__main__':main()
