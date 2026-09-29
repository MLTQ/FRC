"""Conservative bulk flux with separate full Poynting/current transport."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_reconstruction import radial_slopes,axial_increments,radial_pairs
from .mhd_positivity import positive_scale
from .toroidal_flux import fast_speed


def scalar_faces(mesh,value,order,*,upper=None):
    """Positive scalar traces; optional upper bound also contracts only slopes."""
    center=.25*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:,None]
    slope=radial_slopes(value,center) if order==2 else np.zeros_like(value)
    axial=axial_increments(value) if order==2 else np.zeros_like(value)
    lo=slope*(.5*mesh.rf[:-1,None]**2-center);hi=slope*(.5*mesh.rf[1:,None]**2-center)
    minimum=value+np.minimum(lo,hi)-.5*abs(axial)
    factor=positive_scale(value,minimum)
    if upper is not None:
        maximum=value+np.maximum(lo,hi)+.5*abs(axial)
        factor=np.minimum(factor,positive_scale(upper-value,upper-maximum))
    low,high=value+factor*lo,value+factor*hi
    left,right=radial_pairs(low,high)
    low,high=value-.5*factor*axial,value+.5*factor*axial
    return left,right,np.concatenate((high[:,-1:],high),axis=1),np.concatenate((low,low[:,:1]),axis=1)


def physical(q,fraction,current_energy,radius,axis,parameters):
    rho,vr,vz,p,br,bz,omega,g=np.moveaxis(q,-1,0);vt=radius*omega;bt=radius*g
    pe=fraction*p;pi=(1-fraction)*p;internal=pe/(parameters.gamma_e-1)+pi/(parameters.gamma_i-1)
    kinetic=.5*rho*(vr*vr+vz*vz+vt*vt);magnetic=(br*br+bz*bz+bt*bt)/(2*MU0)
    vn,bn=(vr,br) if axis==0 else (vz,bz);angular=rho*radius*radius*omega
    energy=internal+kinetic+magnetic+current_energy
    conserved=np.stack((rho,rho*vr,rho*vz,energy,angular),axis=-1)
    flux=np.stack((rho*vn,rho*vr*vn-br*bn/MU0,rho*vz*vn-bz*bn/MU0,
                   (internal+p+kinetic+current_energy)*vn,angular*vn-radius*bt*bn/MU0),axis=-1)
    flux[...,1+axis]+=p+magnetic
    gamma=parameters.gamma_e*fraction+parameters.gamma_i*(1-fraction)
    speed=abs(vn)+fast_speed(q,radius,axis,gamma)
    return conserved,flux,speed


def rusanov(left,right,fl,fr,kl,kr,normal_field,radius,axis,parameters):
    left,right=left.copy(),right.copy();left[...,4+axis]=right[...,4+axis]=normal_field
    ul,vl,al=physical(left,fl,kl,radius,axis,parameters)
    ur,vr,ar=physical(right,fr,kr,radius,axis,parameters)
    return .5*(vl+vr)-.5*np.maximum(al,ar)[...,None]*(ur-ul)
