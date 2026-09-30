"""Linear radial acoustic spectrum and independent full-RHS perturbation controls."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig,expm,solve
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.log_fluid import LogFluid
from .log_fluid_reference import LogFluidReference
from .log_fluid_records import provenance,sha


def calculate(spec):
    nr=spec['nr'];ref=LogFluidReference(gamma_i=5/3);mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    g=CompatibleGeometry(mesh,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    fluid=LogFluid(g,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    f0=fluid.transport.rfaces[0].toarray();f1=fluid.transport.rfaces[1].toarray();f1[[0,-1]]=0.
    D=np.diff(mesh.rf[:,None]*f1,axis=0)/mesh.radial_measure[:,None]
    R0=fluid.transport.radial[0].toarray();R1=fluid.transport.radial[1].toarray()
    _,weights=np.polynomial.legendre.leggauss(g.order)
    hoop=np.sum(.5*mesh.dr*weights[None,:,None]*R0.reshape(nr,g.order,nr),axis=1)/mesh.radial_measure[:,None]
    G=-np.diff(mesh.rf[:,None]*f0,axis=0)/mesh.radial_measure[:,None]+hoop
    Hn=R0.T@(g.wr[:,None]*R0);Hm=R1.T@(g.wr[:,None]*R1)
    adjoint=solve(Hm,D.T@Hn,assume_a='pos')
    zero=np.zeros_like(D);A=ref.sound*np.block([[zero,-D],[G,zero]])
    stable=ref.sound*np.block([[zero,-D],[adjoint,zero]])
    values,vectors=eig(A);which=np.argmax(values.real);mode=vectors[:,which];mode/=np.max(abs(mode))
    values_stable=eig(stable,right=False);x=mode.real
    if np.linalg.norm(x)<.1:x=mode.imag
    H=np.block([[Hn,zero],[zero,Hm]])
    unweighted_defect=G-adjoint
    defect=np.linalg.norm(H@A+A.T@H)/np.linalg.norm(H@A)
    # Compare the analytical block linearization with central differences of the
    # original nonlinear equations. No eigensolver is used by that full RHS.
    n0=ref.number_density;shape=(nr,mesh.nz);node=(nr+1,mesh.nz)
    entropy=[np.full(node,np.log(p/(n0*E*(n0/1e20)**(gamma-1))))
             for p,gamma in ((.6*ref.pressure,ref.gamma_e),(.4*ref.pressure,ref.gamma_i))]
    def action(delta):
        density=n0*np.broadcast_to((1+delta*x[:nr])[:,None],shape).copy()
        momentum=np.zeros(shape+(3,));momentum[...,0]=n0*ref.sound*delta*x[nr:,None]
        state=fluid.initialize(density,momentum,*entropy);rate,_=fluid.derivative(state)
        return np.r_[rate.density.mean(axis=1)/n0,rate.momentum[...,0].mean(axis=1)/(n0*ref.sound)]
    controls=[]
    for epsilon in (1e-4,1e-5,1e-6):
        actual=(action(epsilon)-action(-epsilon))/(2*epsilon);wanted=A@x
        controls.append(dict(epsilon=epsilon,relative_L2=float(np.linalg.norm(actual-wanted)/np.linalg.norm(wanted))))
    # Exact matrix exponential verifies growth independently of an RK4 step.
    end=3/values[which].real
    exact=expm(A*end)@x;energy0=float(x@H@x);energy=float(exact@H@exact)
    temporal=[]
    for count in (64,128,256):
        state=x.copy();dt=end/count
        for _ in range(count):
            k1=A@state;k2=A@(state+.5*dt*k1);k3=A@(state+.5*dt*k2);k4=A@(state+dt*k3)
            state+=dt*(k1+2*k2+2*k3+k4)/6
        temporal.append(dict(steps=count,relative_error=float(np.linalg.norm(state-exact)/np.linalg.norm(exact)),energy_ratio=float(state@H@state)/energy0))
    arrays=dict(divergence=D,pressure=G,adjoint_control=adjoint,density_metric=Hn,momentum_metric=Hm,
                operator=A,eigenvalues=values,adjoint_eigenvalues=values_stable,fastest_mode=mode,
                initial=x,exact_final=exact)
    record=dict(spec=spec,sound=ref.sound,maximum_growth_per_s=float(values[which].real),
        growth_times_radius_over_sound=float(values[which].real*ref.radius/ref.sound),
        adjoint_maximum_real_per_s=float(max(values_stable.real)),physical_energy_adjoint_defect=float(defect),
        gradient_relative_difference=float(np.linalg.norm(unweighted_defect)/np.linalg.norm(G)),
        nonlinear_difference_controls=controls,end_time=end,exact_energy_ratio=energy/energy0,rk4_controls=temporal,
        density_metric_condition=float(np.linalg.cond(Hn)),momentum_metric_condition=float(np.linalg.cond(Hm)))
    return record,arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);records=[];identity=provenance();probe=Path(__file__)
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32)]+[dict(nr=16,degree=5,quadrature_order=16)]
    for i,spec in enumerate(specs):
        record,arrays=calculate(spec);path=d/f'radial_{i:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record);print(json.dumps(record),flush=True)
    (d/'probe_source.py').write_bytes(probe.read_bytes());(d/'probe_source.md').write_bytes(probe.with_suffix('.md').read_bytes())
    (d/'records.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],probe_sha256=sha(probe),records=records),indent=2)+'\n')


if __name__=='__main__':main()
