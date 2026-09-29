"""Independent parent-species forcing and conservative moving-state reference."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from plasma_bench.finite_mass_state import FiniteMassParameters,from_primitive
from plasma_bench.finite_mass_rhs import FiniteMassSources
from plasma_bench.toroidal_state import radius_squared
from .moving_species_reference import parent_sources


def point_conserved(ref,r,z,t,*,rate=False):
    f=ref.fields(r,z,t);total=ref.ion_mass+ref.electron_mass;rho=total*f['n']
    if not rate:
        energy=f['pe']/(ref.gamma_e-1)+f['pi']/(ref.gamma_i-1)+np.sum(f['b']**2,axis=-1)/(2*MU0)
        for key,mass in (('i',ref.ion_mass),('e',ref.electron_mass)):
            energy+=.5*mass*f['n']*np.sum(f['v'+key]**2,axis=-1)
        return np.stack((rho,rho*f['v'][...,0],rho*f['v'][...,2],energy,rho*r*f['v'][...,1]),axis=-1)
    dt=ref.derivatives(r,z,t)[2];drho=total*dt['n']
    momentum=drho[...,None]*f['v']+rho[...,None]*dt['v']
    energy=dt['pe']/(ref.gamma_e-1)+dt['pi']/(ref.gamma_i-1)+np.sum(f['b']*dt['b'],axis=-1)/MU0
    for key,mass in (('i',ref.ion_mass),('e',ref.electron_mass)):
        energy+=.5*mass*dt['n']*np.sum(f['v'+key]**2,axis=-1)+mass*f['n']*np.sum(f['v'+key]*dt['v'+key],axis=-1)
    return np.stack((drho,momentum[...,0],momentum[...,2],energy,r*momentum[...,1]),axis=-1)


def cell_average(mesh,function,order=8):
    out=None;nodes,weights=np.polynomial.legendre.leggauss(order)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for z,wz in zip(nodes,weights):
            zz=mesh.z[None]+.5*mesh.dz*z;value=function(*np.broadcast_arrays(r,zz))
            weight=.25*wx*wz*r*mesh.dr/mesh.radial_measure[:,None]
            result=value*weight[(...,)+(None,)*(value.ndim-2)]
            out=result if out is None else out+result
    return out


def initialize(mesh,ref,time=0.):
    par=FiniteMassParameters(ref.ion_mass,model='hall' if ref.electron_mass==0 else 'two_fluid',
        gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    exact=cell_average(mesh,lambda r,z:point_conserved(ref,r,z,time))
    pi=cell_average(mesh,lambda r,z:ref.fields(r,z,time)['pi'])
    g=cell_average(mesh,lambda r,z:ref.fields(r,z,time)['g'])
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');pe=ref.fields(r,z,time)['pe']
    psi=ref.fields(*mesh.corners,time)['psi'];psi[:,-1]=psi[:,0]
    rho=exact[...,0]
    state=from_primitive(mesh,par,psi,g,rho,pe,pi,vr=exact[...,1]/rho,vz=exact[...,2]/rho,
        omega=exact[...,4]/(rho*radius_squared(mesh)))
    return state,par


class ManufacturedSources:
    def __init__(self,mesh,reference):self.mesh,self.reference=mesh,reference

    def __call__(self,time,evaluation):
        ref=self.reference;model=evaluation.model;mesh=self.mesh
        def fluid(r,z):
            forces,heat,_=parent_sources(ref,r,z,time);fields=ref.fields(r,z,time)
            force=forces['i']+forces['e'];power=heat['i']+heat['e']
            power+=np.sum(fields['vi']*forces['i']+fields['ve']*forces['e'],axis=-1)
            return np.stack((np.zeros_like(r),force[...,0],force[...,2],power,r*force[...,1]),axis=-1)
        source=cell_average(mesh,fluid)
        r,z=model.geometry.r,model.geometry.z;forces,_,_=parent_sources(ref,r,z,time)
        density=ref.fields(r,z,time)['n'];mass=ref.ion_mass+ref.electron_mass
        electric=(ref.ion_mass*forces['e']-ref.electron_mass*forces['i'])/(mass*E*density[...,None])
        # Axisymmetric scalar heat has a regular limit; evaluate r=0 with the
        # analytic even extension sampled by complex-step-safe limiting r.
        r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');r=r.copy();r[0]=1e-10*ref.radius
        heat=parent_sources(ref,r,z,time)[1]['e'];n=ref.fields(r,z,time)['n']
        entropy=(ref.gamma_e-1)*heat/(n/model.electron.density_scale)**(ref.gamma_e-1)
        return FiniteMassSources(source,entropy,electric)
