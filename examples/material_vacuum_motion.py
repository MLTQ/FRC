"""Evolve the magnetic exterior along a prescribed nonaffine boundary motion."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_magnetic_vacuum import MaterialMagneticVacuum
from .implicit_plasma_records import provenance,freeze,sha


def dipole(r,z,moment):
    square=r*r+z*z;scale=moment/square**2.5
    return moment*r*r/square**1.5,np.column_stack((3*scale*r*z,scale*(2*z*z-r*r)))


def reference(coordinates,velocity,ids,moment):
    """Independent dipole Green identity, Maxwell traction and domain energy rate."""
    q,w=np.polynomial.legendre.leggauss(64);u=.5*(q+1);w=.5*w
    boundary=np.r_[ids[:,0],ids[-1,1:],ids[-2::-1,-1]]
    gradient=np.zeros_like(coordinates);energy=0.;power=0.
    for a,b in zip(boundary[:-1],boundary[1:]):
        tangent=coordinates[b]-coordinates[a];normal=np.array([tangent[1],-tangent[0]])
        x=(1-u[:,None])*coordinates[a]+u[:,None]*coordinates[b]
        v=(1-u[:,None])*velocity[a]+u[:,None]*velocity[b]
        psi,field=dipole(x[:,0],x[:,1],moment);square=np.sum(field*field,axis=1)
        energy-=np.pi/MU0*np.sum(w*psi*(field[:,1]*normal[0]-field[:,0]*normal[1]))
        load=2*np.pi/MU0*(w*x[:,0])[:,None]*(.5*square[:,None]*normal-(field@normal)[:,None]*field)
        gradient[a]+=np.sum((1-u[:,None])*load,axis=0);gradient[b]+=np.sum(u[:,None]*load,axis=0)
        power-=np.pi/MU0*np.sum(w*x[:,0]*square*(v@normal))
    gradient[coordinates[:,0]==0,0]=0
    return dict(energy=float(energy),gradient=gradient,power=float(power))


def run(destination,*,nr=4,degree=3,quadrature=6,vacuum_cells=4,modes=6,steps=20):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh physical trajectory directory required')
    mesh=MHDMesh(.035,.07,nr,2*nr);ids=np.arange((nr+1)*(2*nr+1)).reshape(nr+1,2*nr+1)
    started=time.monotonic();identity=provenance(['examples.material_vacuum_motion']);out.mkdir(parents=True);freeze(out,identity)
    record=dict(terminal=False,success=False,completed_steps=0,steps=steps,nr=nr,nz=2*nr,degree=degree,quadrature=quadrature,
        vacuum_cells=vacuum_cells,modes=modes,source_sha256=identity['source_sha256'],duration=7.121769264796269e-6,records=[],
        scope='Prescribed nonaffine magnetic-only exterior motion with an independent fixed physical dipole. No evolved fluid, coil-driven formation or laboratory-validity claim.')
    def save():
        record.update(elapsed_seconds=time.monotonic()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        model=MaterialMagneticVacuum(mesh,degree=degree,quadrature_order=quadrature,vacuum_cells=vacuum_cells,modes=modes)
        r,z=model.reference.T;displacement=np.column_stack((.25*r*(1-(z/mesh.half_length)**2),.15*z*(1-(r/mesh.radius)**2)))
        velocity=displacement/record['duration'];moment=2e-4*mesh.radius**3
        for step in range(steps+1):
            coordinates=model.reference+(step/steps)*displacement
            x=coordinates[model.boundary_nodes];boundary,field=dipole(x[:,0],x[:,1],moment)
            flux_rate=np.sum(np.column_stack((x[:,0]*field[:,1],-x[:,0]*field[:,0]))*velocity[model.boundary_nodes],axis=1)
            actual=model.evaluate(coordinates,boundary,velocity);exact=reference(coordinates,velocity,ids,moment)
            points=actual['points'];_,target=dipole(points[:,0],points[:,1],moment);volume=actual['volume']
            rate=actual['geometry_power']+actual['boundary_load']@flux_rate
            row=dict(step=step,time=step/steps*record['duration'],energy_j=actual['energy'],reference_energy_j=exact['energy'],
                energy_relative_error=actual['energy']/exact['energy']-1,
                field_relative_l2_error=float(np.sqrt(np.sum(volume[:,None]*(actual['magnetic_field']-target)**2)/np.sum(volume[:,None]*target**2))),
                traction_nodal_relative_l2_error=float(np.linalg.norm(actual['coordinate_energy_gradient']-exact['gradient'])/np.linalg.norm(exact['gradient'])),
                geometry_power_w=actual['geometry_power'],matrix_geometry_power_w=float(.5*boundary@actual['hessian_rate']@boundary),
                total_power_w=float(rate),reference_power_w=exact['power'],power_relative_error=float(rate/exact['power']-1),
                minimum_vacuum_map_jacobian=float(actual['determinant'].min()))
            filename=out/f'step_{step:05}.npz'
            np.savez_compressed(filename,coordinates=coordinates,boundary_flux=boundary,points=points,field=actual['magnetic_field'],
                reference_field=target,force=actual['force'],reference_force=-exact['gradient'],hessian=actual['hessian'],hessian_rate=actual['hessian_rate'])
            row.update(arrays_file=filename.name,arrays_sha256=sha(filename));record['records'].append(row);record['completed_steps']=step;save()
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3 or shutil.disk_usage('.').free<256*1024**2:raise RuntimeError('Own4GiB/256MiB resource guard')
            print('material magnetic vacuum',nr,step,'/',steps,flush=True)
        record['success']=True
    except Exception as exc:record.update(failure_type=type(exc).__name__,failure_message=str(exc))
    record['terminal']=True
    if provenance(['examples.material_vacuum_motion'])!=identity:raise RuntimeError('Source changed during run')
    save();print('material magnetic vacuum terminal',record['success'],record.get('failure_message',''),flush=True)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--nr',type=int,default=4);p.add_argument('--degree',type=int,default=3)
    p.add_argument('--quadrature',type=int,default=6);p.add_argument('--vacuum-cells',type=int,default=4);p.add_argument('--modes',type=int,default=6);p.add_argument('--steps',type=int,default=20)
    a=p.parse_args();d=run(a.destination,nr=a.nr,degree=a.degree,quadrature=a.quadrature,vacuum_cells=a.vacuum_cells,modes=a.modes,steps=a.steps)
    raise SystemExit(0 if d['success'] else 1)
