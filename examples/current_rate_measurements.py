"""Physical instantaneous Bdot after inertia, with an independent spectral oracle."""
from pathlib import Path
from dataclasses import asdict
import hashlib
import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import splu
from scipy.constants import atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from plasma_bench.hall_matrices import transport_matrix
from .electron_reference import ElectronReference
from .current_inertia_reference import g_rate
from .current_rate_reference import CurrentRateReference
from .current_thermal_records import provenance
from .verification_records import write_record


def discrete(model,psi,g,ue,ui):
    geo=model.geometry;solver=CurrentThermalMidpoint(model);packed=geo.pack(psi,g);j=model.current(packed)
    star=model.star_field(packed,j);correction=star-geo.magnetic_samples(psi,g)
    pressure=model.pressure(ue,ui)[0]
    full=model.delta*transport_matrix(geo,star,model.density)@j+model.mass@pressure
    conv=model.delta*transport_matrix(geo,correction,model.density)@j
    # tau is only a variable scaling for the static linear equation, not a step.
    tau=1e-8;impedance=float(model.inertia_coefficient.max())/tau;d=solver.diagonal
    matrix=bmat([[solver.scaled_mass,tau*impedance*d@solver.stiffness@d],
                 [-d@model.inertia@d/(tau*impedance),solver.scaled_mass]],format='csc')
    rhs=np.column_stack([np.r_[np.zeros_like(j),solver.scale*f/impedance] for f in (full,conv)])
    solved=splu(matrix).solve(rhs);residual=np.linalg.norm(matrix@solved-rhs,axis=0)/np.linalg.norm(rhs,axis=0)
    if residual.max()>1e-10:raise RuntimeError('Static current-inertia solve did not meet residual')
    rates=[]
    for col in range(2):
        electric=impedance*solver.scale*solved[len(j):,col]
        rates.append(geo.unpack(geo.faraday@electric)[1])
    return rates,residual.tolist()


def run(directory=Path('outputs/current_rate_verification'),method='mixed'):
    if method not in ('mixed','schur'):raise ValueError('Choose mixed or schur rate solve')
    solver=discrete
    if method=='schur':
        from .current_rate_schur import discrete as solver
    directory.mkdir(parents=True,exist_ok=True);ref=ElectronReference();reference_cases=[];references=[]
    names=('examples/current_rate_reference.py','examples/current_rate_measurements.py','examples/current_inertia_reference.py')
    if method=='schur':names+=('examples/current_rate_schur.py',)
    source=dict(engine=provenance(),rate_sources={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names})
    for order in (64,128,256):
        oracle=CurrentRateReference(ref,order,256)
        rr=ref.radius*np.sqrt(oracle.s[:,None]);rr[0]=ref.radius*1e-8
        conv=g_rate(ref,rr,oracle.z[None],oracle.mu,oracle.delta)
        values=[];diagnostics=[]
        for forcing in (oracle.original_source(),conv):
            value,diagnostic=oracle.solve(forcing);values.append(value);diagnostics.append(diagnostic)
        probe=oracle.original_source(axis_probe_fraction=1e-6)
        probe_error=float(np.linalg.norm(probe-oracle.original_source())/np.linalg.norm(probe))
        path=directory/f'continuum_{order}.npz';np.savez_compressed(path,s=oracle.s,z=oracle.z,full=values[0],convective=values[1])
        reference_cases.append(dict(radial_order=order,axial_points=256,solves=diagnostics,axis_probe_relative_change=probe_error,
            fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())))
        references.append((oracle,values));print('continuum',order,diagnostics,flush=True)
    records=[]
    for nr in (8,16,32,64,128):
        mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr);f=ref.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))
        model=CurrentThermal(mesh,f['density'],39.948*atomic_mass);psi=np.zeros((nr+1,2*nr+1));g=ref.cell_average(mesh,'g')
        ue=model.electron.from_pressure(f['pressure']);ui=model.ion.from_pressure(np.full_like(ue,16.))
        values,residual=solver(model,psi,g,ue,ui);metrics={};arrays=dict(radial_faces=mesh.rf,axial_faces=mesh.zf)
        for i,key in enumerate(('full','convective')):
            targets=[oracle.cell_average(mesh,v[i]) for oracle,v in references];target=targets[-1];error=values[i]-target
            metrics[key]=dict(relative_L2=float(np.sqrt(np.sum(mesh.volumes*error**2)/np.sum(mesh.volumes*target**2))),
                axis_relative_L2=float(np.linalg.norm(error[0])/np.linalg.norm(target[0])),
                maximum_absolute_error=float(abs(error).max()),maximum_reference=float(abs(target).max()),
                continuum_64_to_128_relative=float(np.linalg.norm(targets[0]-targets[1])/np.linalg.norm(targets[1])),
                continuum_128_to_256_relative=float(np.linalg.norm(targets[1]-target)/np.linalg.norm(target)))
            arrays.update({key+'_computed':values[i],key+'_reference':target})
        path=directory/f'rate_{nr}.npz';np.savez_compressed(path,**arrays)
        records.append(dict(nr=nr,nz=2*nr,linear_residuals=residual,metrics=metrics,
            fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())))
        print('rate',nr,metrics,flush=True)
    assert source['engine']==provenance()
    assert source['rate_sources']=={n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}
    write_record(directory/'measurements.json',dict(completed=True,method=method,provenance=source,reference=asdict(ref),
        ion_mass=39.948*atomic_mass,ion_pressure=16.,continuum=reference_cases,cases=records,
        scope='Independent physical instantaneous gdot after inertia inversion; full forcing and isolated convective contribution are distinct.'))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--method',choices=('mixed','schur'),default='mixed')
    parser.add_argument('--directory',default='outputs/current_rate_verification');args=parser.parse_args()
    run(Path(args.directory),args.method)
