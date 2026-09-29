"""Independent volume initialization and rate/evolution norms for coupled waves."""
import numpy as np
from plasma_bench.finite_mass_state import FiniteMassParameters,from_primitive,evaluate
from plasma_bench.toroidal_state import radius_squared
from .finite_mass_manufactured import cell_average
from .current_thermal_norms import vectors,scalar


def initialize(mesh,ref):
    par=FiniteMassParameters(ref.ion_mass,model='hall' if ref.electron_mass==0 else 'two_fluid',gamma_e=ref.gamma,gamma_i=ref.gamma_i)
    def matter(r,z):
        f=ref.fields(r,z);rho=ref.density+f['density'];v=f['velocity']
        return np.stack((rho,rho*v[...,0],rho*v[...,2],rho*r*v[...,1]),axis=-1)
    matter=cell_average(mesh,matter);g=cell_average(mesh,lambda r,z:ref.fields(r,z)['magnetic'][...,1]/r)
    r,z=mesh.corners;psi=.5*ref.bias*r*r+ref.fields(r,z)['psi'];psi[:,-1]=psi[:,0]
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');pe=ref.electron_pressure+ref.fields(r,z)['electron_pressure']
    pi=ref.ion_pressure+cell_average(mesh,lambda r,z:ref.fields(r,z)['ion_pressure'])
    rho=matter[...,0]
    return from_primitive(mesh,par,psi,g,rho,pe,pi,vr=matter[...,1]/rho,vz=matter[...,2]/rho,
        omega=matter[...,3]/(rho*radius_squared(mesh))),par


def norms(mesh,ref,time,state,parameters,*,order=8):
    ev=evaluate(mesh,state,parameters);q=ev.primitive;model=ev.model
    pe=model.electron.pressure(state.electron_entropy);pi=(parameters.gamma_i-1)*ev.components['ion']
    errors={k:0. for k in ('magnetic','current','density','velocity','electron_pressure','ion_pressure')};scales=errors.copy()
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for z,wz in zip(nodes,weights):
            az=.5*(z+1);zz=mesh.z[None]+.5*mesh.dz*z;f=ref.fields(r,zz,time)
            b,j,ax=vectors(mesh,state.psi,state.g,ev.current,r,az);b[...,2]-=ref.bias
            v=np.stack(np.broadcast_arrays(q[...,1],r*q[...,6],q[...,2]),axis=-1)
            values=dict(magnetic=b,current=j,density=state.conserved[...,0]-ref.density,velocity=v,
                electron_pressure=scalar(pe,ax,az)-ref.electron_pressure,ion_pressure=pi-ref.ion_pressure)
            weight=.5*np.pi*wx*wz*r*mesh.dr*mesh.dz
            for key,value in values.items():
                error=(value-f[key])**2;norm=f[key]**2
                if value.ndim==3:error=error.sum(axis=-1);norm=norm.sum(axis=-1)
                errors[key]+=np.sum(weight*error);scales[key]+=np.sum(weight*norm)
    return {key:dict(relative_L2=float(np.sqrt(errors[key]/scales[key])),rms_error=float(np.sqrt(errors[key]/np.sum(mesh.volumes)))) for key in errors}
