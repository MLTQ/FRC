"""Nonlinear inertia forcing convergence, separate from the Beltrami cancellation."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from .electron_reference import ElectronReference
from .current_inertia_reference import fields,g_rate
from .current_thermal_records import provenance
from .verification_records import write_record


def measure(nr):
    ref=ElectronReference();mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr)
    n=ref.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))['density']
    model=CurrentThermal(mesh,n,39.948*atomic_mass);geo=model.geometry
    psi=np.zeros((nr+1,2*nr+1));g=ref.cell_average(mesh,'g');packed=geo.pack(psi,g);j=model.current(packed)
    correction=model.star_field(packed,j)-geo.magnetic_samples(psi,g)
    electric=model.delta*np.cross(geo.sample_electric(j)/(model.density[...,None]*E),correction)
    coefficients=geo.project_electric(electric);rate=geo.unpack(geo.faraday@coefficients)[1]
    exact=fields(ref,geo.r,geo.z,model.reduced_mass,model.delta)
    nodes,weights=np.polynomial.legendre.leggauss(12);target=np.zeros_like(g)
    lo,hi=.5*mesh.rf[:-1]**2,.5*mesh.rf[1:]**2
    for a,wa in zip(nodes,weights):
        r=np.sqrt(lo+hi+a*(hi-lo))[:,None]
        for b,wb in zip(nodes,weights):
            z=mesh.z[None]+.5*mesh.dz*b
            target+=.25*wa*wb*g_rate(ref,r,z,model.reduced_mass,model.delta)
    def norm(error,scale,weight):return float(np.sqrt(np.sum(weight*error*error)/np.sum(weight*scale*scale)))
    corr=exact['star_correction'];ee=exact['reduced_convective_electric']
    metrics=dict(star_correction_relative_L2=norm(correction-corr,corr,geo.weight[...,None]),
        axis_star_correction_relative_L2=norm(correction[0]-corr[0],corr[0],geo.weight[0,...,None]),
        electric_relative_L2=norm(electric-ee,ee,geo.weight[...,None]),
        g_forcing_relative_L2=norm(rate-target,target,mesh.volumes),
        axis_g_forcing_relative_L2=norm(rate[0]-target[0],target[0],mesh.volumes[0]),
        g_forcing_max_absolute_error=float(abs(rate-target).max()),g_forcing_max_reference=float(abs(target).max()))
    return dict(nr=nr,nz=2*nr,reference=asdict(ref),ion_mass=model.ion_mass,electron_mass=model.electron_mass,
        reduced_mass=model.reduced_mass,delta=model.delta,quadrature=5,reference_quadrature=12,metrics=metrics),dict(
        g=g,number_density=n,computed_g_forcing=rate,reference_g_forcing=target,radial_faces=mesh.rf,axial_faces=mesh.zf)


def run(directory=Path('outputs/current_inertia_spatial')):
    directory.mkdir(parents=True,exist_ok=True)
    names=('examples/current_inertia_reference.py','examples/current_inertia_spatial.py')
    source=dict(engine=provenance(),spatial_sources={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names})
    records=[]
    for nr in (8,16,32,64,128):
        record,arrays=measure(nr);path=directory/f'nonlinear_{nr}.npz';np.savez_compressed(path,**arrays)
        record.update(fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        records.append(record);print(nr,record['metrics'],flush=True)
    assert source['engine']==provenance()
    assert source['spatial_sources']=={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}
    write_record(directory/'measurements.json',dict(completed=True,provenance=source,cases=records,
        interpretation='Instantaneous nonlinear convective forcing of B*, not Bdot after the inertial solve and not an evolution.'))


if __name__=='__main__':run()
