"""Independent initialization and physical norms for moving-current rate cases."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from .moving_species_reference import MovingSpeciesReference,parent_sources
from .current_thermal_norms import vectors


def initialize(case):
    ref=MovingSpeciesReference(**case['reference']);t=case['phase']/ref.frequency
    mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],2*case['nr'])
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=ref.fields(r,z,t)
    model=CurrentThermal(mesh,f['n'],ref.ion_mass,model='hall' if ref.electron_mass==0 else 'two_fluid',
        gamma_e=ref.gamma_e,gamma_i=ref.gamma_i,quadrature_order=case['quadrature'])
    psi=ref.fields(*mesh.corners,t)['psi'];psi[:,-1]=psi[:,0]
    g=np.zeros((mesh.nr,mesh.nz));nodes,weights=np.polynomial.legendre.leggauss(12)
    for x,wx in zip(nodes,weights):
        rr=mesh.r[:,None]+.5*mesh.dr*x
        for q,wz in zip(nodes,weights):
            g+=.25*wx*wz*rr*mesh.dr/mesh.radial_measure[:,None]*ref.fields(rr,mesh.z[None]+.5*mesh.dz*q,t)['g']
    ue=model.electron.from_pressure(f['pe']);ui=model.ion.from_pressure(f['pi']);nt=ref.derivatives(r,z,t)[2]['n']
    geo=model.geometry;force=parent_sources(ref,geo.r,geo.z,t)[0];physical=ref.fields(geo.r,geo.z,t)
    source=(model.beta*force['e']-model.alpha*force['i'])/(E*physical['n'][...,None])
    return ref,model,psi,g,ue,ui,f['regular_velocity'],nt,source,t


def norms(mesh,ref,time,rate_psi,rate_g,rate_current,electric,*,order=8):
    errors={k:0. for k in ('magnetic_rate','current_rate','electric')};scales=errors.copy()
    axis=errors.copy();axis_scale=errors.copy();maximum=errors.copy();reference_max=errors.copy()
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for q,wz in zip(nodes,weights):
            z=mesh.z[None]+.5*mesh.dz*q;a=ref.fields(r,z,time);dt=ref.derivatives(r,z,time)[2]
            b,j,_=vectors(mesh,rate_psi,rate_g,rate_current,r,.5*(1+q))
            _,ef,_=vectors(mesh,rate_psi,rate_g,electric,r,.5*(1+q))
            weight=.5*np.pi*wx*wz*r*mesh.dr*mesh.dz
            for key,value,target in (('magnetic_rate',b,dt['b']),('current_rate',j,dt['j']),('electric',ef,a['e'])):
                error=np.sum((value-target)**2,axis=-1);square=np.sum(target**2,axis=-1)
                errors[key]+=np.sum(weight*error);scales[key]+=np.sum(weight*square)
                axis[key]+=np.sum(weight[0]*error[0]);axis_scale[key]+=np.sum(weight[0]*square[0])
                maximum[key]=max(maximum[key],float(np.sqrt(error.max())))
                reference_max[key]=max(reference_max[key],float(np.sqrt(square.max())))
    return {key:dict(relative_L2=float(np.sqrt(errors[key]/scales[key])),
        axis_relative_L2=float(np.sqrt(axis[key]/axis_scale[key])),
        rms_error=float(np.sqrt(errors[key]/np.sum(mesh.volumes))),
        maximum_error=maximum[key],maximum_reference=reference_max[key]) for key in errors}
