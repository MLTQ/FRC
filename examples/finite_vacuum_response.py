"""Independent finite-cylinder exterior dipole comparison, without plasma dynamics."""
import argparse
import json
from pathlib import Path
import time
import resource
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_cylinder_vacuum import FiniteCylinderVacuum
from .implicit_plasma_records import provenance,freeze,sha


def dipole(r,z,strength=1e-8,offset=.1):
    r,z=np.broadcast_arrays(np.asarray(r,float),np.asarray(z,float)-offset);rr=np.hypot(r,z)
    return strength*r*r/rr**3,3*strength*r*z/rr**5,strength*(2*z*z-r*r)/rr**5


def run(destination, *, nr=8,nz=16,vacuum_cells=8,modes=8):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh result directory required')
    identity=provenance(['examples.finite_vacuum_response']);out.mkdir(parents=True);freeze(out,identity);started=time.monotonic()
    mesh=MHDMesh(.075,.5,nr,nz)
    vacuum=FiniteCylinderVacuum(mesh,vacuum_cells=vacuum_cells,modes=modes)
    boundary=dipole(vacuum.boundary_r,vacuum.boundary_z)[0]
    # Independent Green boundary integral on the inner cylinder; normals point
    # from plasma into vacuum, opposite to the exterior volume's own normal.
    q,w=np.polynomial.legendre.leggauss(160);energy=0.
    for r,z,weights,normal in ((np.full(160,mesh.radius),mesh.half_length*q,mesh.half_length*w,0),
        (.5*mesh.radius*(q+1),np.full(160,mesh.half_length),.5*mesh.radius*w,1),
        (.5*mesh.radius*(q+1),np.full(160,-mesh.half_length),.5*mesh.radius*w,-1)):
        psi,br,bz=dipole(r,z);normal_derivative=r*bz if normal==0 else -normal*r*br
        energy-=np.pi/MU0*np.sum(weights/r*psi*normal_derivative)
    rr=np.geomspace(mesh.radius,3*mesh.radius,25);zz=np.linspace(-mesh.half_length,mesh.half_length,101)
    r,z=np.meshgrid(rr,zz,indexing='ij');actual=vacuum.samples(boundary,r,z);expected=dipole(r,z)
    err=np.sqrt(np.mean((actual[1]-expected[1])**2+(actual[2]-expected[2])**2))
    scale=np.sqrt(np.mean(expected[1]**2+expected[2]**2))
    # Real end-cap/axis points and points beyond the multipole matching box.
    end_r=np.linspace(0,mesh.radius,31);end_z=np.array([-.65,-.5,.5,.65,1.2])
    er,ez=np.meshgrid(end_r,end_z,indexing='ij');end=vacuum.samples(boundary,er,ez);target=dipole(er,ez)
    end_error=np.sqrt(np.mean((end[1]-target[1])**2+(end[2]-target[2])**2))/np.sqrt(np.mean(target[1]**2+target[2]**2))
    coefficients=vacuum.coefficients(boundary);native=float(.5*boundary@vacuum.hessian@boundary)
    np.savez_compressed(out/'field.npz',boundary=boundary,coefficients=coefficients,r=r,z=z,psi=actual[0],br=actual[1],bz=actual[2],
        exact_br=expected[1],exact_bz=expected[2],end_r=er,end_z=ez,end_br=end[1],end_bz=end[2])
    report=dict(terminal=True,nr=nr,nz=nz,vacuum_cells=vacuum_cells,modes=modes,
        energy_j=native,reference_energy_j=float(energy),energy_relative_difference=float(native/energy-1),
        side_field_relative_rms=float(err/scale),end_field_relative_rms=float(end_error),
        full_energy_j=float(.5*coefficients@vacuum.full_hessian@coefficients),
        outer_radius_m=vacuum.outer.radius,outer_half_length_m=vacuum.outer.half_length,
        source_sha256=identity['source_sha256'],arrays_sha256=sha(out/'field.npz'),elapsed_seconds=time.monotonic()-started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Static magnetic vacuum response to the independently prescribed trace of an off-center axial dipole. No plasma evolution or finite-coil integration yet.')
    if provenance(['examples.finite_vacuum_response'])!=identity:raise RuntimeError('Source changed during finite vacuum solve')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--nr',type=int,default=8);parser.add_argument('--nz',type=int,default=16)
    parser.add_argument('--vacuum-cells',type=int,default=8);parser.add_argument('--modes',type=int,default=8);args=parser.parse_args()
    run(args.destination,nr=args.nr,nz=args.nz,vacuum_cells=args.vacuum_cells,modes=args.modes)
