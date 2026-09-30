"""Independent cylindrical acoustic frequencies and original nonlinear RHS controls."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig,expm,solve
from scipy.special import jn_zeros
from scipy.constants import elementary_charge as E
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.face_radial_fluid import FaceRadialFluid
from .log_fluid_reference import LogFluidReference
from .face_radial_records import provenance,sha


def calculate(spec):
    nr=spec['nr'];ref=LogFluidReference(gamma_i=5/3);mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    g=CompatibleGeometry(mesh,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    fluid=FaceRadialFluid(g,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    R=fluid.transport.radial[0].toarray();Q=ref.radius*fluid.kinematics.radial.toarray()
    D=ref.radius*np.diff(np.pad(np.eye(nr-1),((1,1),(0,0))),axis=0)/mesh.radial_measure[:,None]
    Hn=R.T@(g.wr[:,None]*R);Hq=Q.T@(g.wr[:,None]*Q)
    G=solve(Hq,D.T@Hn,assume_a='pos')
    A=ref.sound*np.block([[np.zeros((nr,nr)),-D],[G,np.zeros((nr-1,nr-1))]])
    H=np.block([[Hn,np.zeros_like(D)],[np.zeros_like(D.T),Hq]])
    values,vectors=eig(A);positive=np.flatnonzero(values.imag>1e-5);positive=positive[np.argsort(values[positive].imag)]
    analytic=ref.sound*jn_zeros(1,3)/ref.radius;measured=values[positive[:3]].imag
    mode=vectors[:,positive[0]];x=mode.real if np.linalg.norm(mode.real)>.1 else mode.imag;x=x/np.max(abs(x))
    n0=ref.number_density;shape=(nr,mesh.nz);node=(nr+1,mesh.nz)
    entropy=[np.full(node,np.log(p/(n0*E*(n0/1e20)**(gamma-1))))
             for p,gamma in ((.6*ref.pressure,ref.gamma_e),(.4*ref.pressure,ref.gamma_i))]
    def action(delta):
        n=n0*np.broadcast_to((1+delta*x[:nr])[:,None],shape).copy()
        q=np.broadcast_to(n0*ref.sound*ref.radius*delta*x[nr:,None],(nr-1,mesh.nz)).copy()
        state=fluid.initialize(n,q,np.zeros(shape+(2,)),*entropy);rate,_=fluid.derivative(state)
        return np.r_[rate.density.mean(axis=1)/n0,rate.radial.mean(axis=1)/(n0*ref.sound*ref.radius)]
    controls=[]
    for epsilon in (1e-4,1e-5,1e-6):
        actual=(action(epsilon)-action(-epsilon))/(2*epsilon);wanted=A@x
        controls.append(dict(epsilon=epsilon,relative_L2=float(np.linalg.norm(actual-wanted)/np.linalg.norm(wanted))))
    end=4*np.pi/measured[0];exact=expm(A*end)@x;energy0=float(x@H@x);temporal=[]
    # Explicit time choices obey the largest eigenfrequency, not only the
    # low-frequency mode selected for comparison.
    base=4*int(np.ceil(end*max(abs(values.imag))/.4/4))
    for count in (base,2*base,4*base):
        y=x.copy();dt=end/count
        for _ in range(count):
            k1=A@y;k2=A@(y+.5*dt*k1);k3=A@(y+.5*dt*k2);k4=A@(y+dt*k3);y+=dt*(k1+2*k2+2*k3+k4)/6
        temporal.append(dict(steps=count,relative_error=float(np.linalg.norm(y-exact)/np.linalg.norm(exact)),energy_ratio=float(y@H@y)/energy0))
    arrays=dict(operator=A,energy_metric=H,eigenvalues=values,initial=x,exact_final=exact,analytic_frequencies=analytic,measured_frequencies=measured)
    record=dict(spec=spec,maximum_real_per_s=float(max(values.real)),maximum_real_scaled=float(max(values.real)*ref.radius/ref.sound),
                physical_energy_adjoint_defect=float(np.linalg.norm(H@A+A.T@H)/np.linalg.norm(H@A)),
                first_three_frequency_relative_errors=((measured-analytic)/analytic).tolist(),
                nonlinear_difference_controls=controls,exact_energy_ratio=float(exact@H@exact)/energy0,rk4_controls=temporal)
    return record,arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);identity=provenance();source=Path(__file__);records=[]
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32)]+[dict(nr=16,degree=5,quadrature_order=16)]
    for index,spec in enumerate(specs):
        record,arrays=calculate(spec);path=d/f'spectrum_{index:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record);print(json.dumps(record),flush=True)
    (d/'probe_source.py').write_bytes(source.read_bytes());(d/'probe_source.md').write_bytes(source.with_suffix('.md').read_bytes())
    (d/'records.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],probe_sha256=sha(source),records=records),indent=2)+'\n')


if __name__=='__main__':main()
