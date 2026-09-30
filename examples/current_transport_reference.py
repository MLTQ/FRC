"""Independent original covariant transport and physical drift-stress divergence."""
import numpy as np
from scipy.constants import elementary_charge as E


def fields(ref,r,z):
    f=ref.fields(r,z);n=f['n'];j=f['current'];y=(r/ref.radius)**2;phase=np.pi*z/ref.half_length
    v=f['particle']/n[...,None];v[...,1]=r*.1*ref.speed/ref.radius*(1+.2*y)*np.cos(phase)
    return n,j,v,j/(E*n[...,None]),ref.reduced_mass/E**2*j[..., :,None]*j[...,None,:]/n[...,None,None]


def exact(ref,r,z,delta):
    r,z=np.broadcast_arrays(r,z);n,j,v,w,stress=fields(ref,r,z);h=1e-24
    dr=[np.imag(x)/h for x in fields(ref,r+1j*h,z)];dz=[np.imag(x)/h for x in fields(ref,r,z+1j*h)]
    def gradient(value,radial,axial):
        result=np.zeros(value.shape+(3,));result[...,0]=radial;result[...,2]=axial
        np.divide(-value[...,1],r,out=result[...,0,1],where=r>0);np.divide(value[...,0],r,out=result[...,1,1],where=r>0)
        result[...,0,1]=np.where(r==0,-radial[...,1],result[...,0,1]);result[...,1,1]=np.where(r==0,radial[...,0],result[...,1,1])
        return result
    vg=gradient(v,dr[2],dz[2]);wg=gradient(w,dr[3],dz[3])
    conn=np.stack((stress[...,0,0]-stress[...,1,1],stress[...,1,0]+stress[...,0,1],stress[...,2,0]),axis=-1)
    np.divide(conn,r[...,None],out=conn,where=r[...,None]>0);conn[...,2]=np.where(r==0,dr[4][...,2,0],conn[...,2])
    force=-(dr[4][..., :,0]+dz[4][..., :,2]+conn)
    return dict(velocity=v,velocity_gradient=vg,drift=w,stress=stress,force=force,
        advection=ref.reduced_mass/E*np.einsum('...ij,...j->...i',wg,v),
        stretching=ref.reduced_mass/E*np.einsum('...ij,...j->...i',vg,w),
        self_drift=-delta*ref.reduced_mass/E*np.einsum('...ij,...j->...i',wg,w))
