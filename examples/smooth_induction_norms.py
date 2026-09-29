"""Independent physical magnetic/current/electric rate norms from saved fields."""
import numpy as np
from plasma_bench.finite_mass_state import evaluate
from .current_thermal_norms import vectors


def measure(mesh,ref,state,par,arrays,time,*,family,order=8):
    ev=evaluate(mesh,state,par);model=ev.model
    current=model.current(model.geometry.pack(arrays['rate_psi'],arrays['rate_g']))
    errors={key:0. for key in ('magnetic_rate','current_rate','electric')};scales=errors.copy()
    axis=errors.copy();axis_scales=errors.copy();maximum=errors.copy()
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for z,wz in zip(nodes,weights):
            az=.5*(1+z);zz=mesh.z[None]+.5*mesh.dz*z
            bd,jd,_=vectors(mesh,arrays['rate_psi'],arrays['rate_g'],current,r,az)
            _,electric,_=vectors(mesh,np.zeros_like(state.psi),np.zeros_like(state.g),arrays['electric'],r,az)
            if family=='manufactured':
                dt=ref.derivatives(r,zz,time)[2];f=ref.fields(r,zz,time)
                exact=dict(magnetic_rate=dt['b'],current_rate=dt['j'],electric=f['e'])
            elif family=='wave':
                dt=ref.derivatives(r,zz,time)[2];f=ref.fields(r,zz,time)
                exact=dict(magnetic_rate=dt['magnetic'].real,current_rate=dt['current'].real,electric=f['electric'])
            else:raise ValueError('Nonzero analytic rate reference required')
            weight=.5*np.pi*wx*wz*r*mesh.dr*mesh.dz
            for key,value in dict(magnetic_rate=bd,current_rate=jd,electric=electric).items():
                squared=np.sum((value-exact[key])**2,axis=-1);scale=np.sum(exact[key]**2,axis=-1)
                errors[key]+=float(np.sum(weight*squared));scales[key]+=float(np.sum(weight*scale))
                axis[key]+=float(np.sum(weight[0]*squared[0]));axis_scales[key]+=float(np.sum(weight[0]*scale[0]))
                maximum[key]=max(maximum[key],float(np.sqrt(squared).max()))
    return {key:dict(relative_L2=float(np.sqrt(errors[key]/scales[key])),
        axis_relative_L2=float(np.sqrt(axis[key]/axis_scales[key])),
        rms_error=float(np.sqrt(errors[key]/mesh.volumes.sum())),maximum_sampled_error=maximum[key]) for key in errors}
