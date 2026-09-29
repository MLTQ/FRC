"""Regular Q1 bulk-velocity sampling and physical cylindrical vorticity."""
import numpy as np


def scalar(geometry,value):
    mesh=geometry.mesh;value=np.asarray(value,float)
    if value.shape!=(mesh.nr+1,mesh.nz) or not np.isfinite(value).all():
        raise ValueError('Finite scalar on complete periodic dual nodes required')
    x=(geometry.r**2-mesh.rf[:-1,None,None,None]**2)/(2*mesh.radial_measure[:,None,None,None])
    z=geometry.az;roll=np.roll(value,-1,axis=1)
    lo=(1-z)*value[:-1,:,None,None]+z*roll[:-1,:,None,None]
    hi=(1-z)*value[1:,:,None,None]+z*roll[1:,:,None,None]
    dx=(hi-lo)/mesh.radial_measure[:,None,None,None]
    dz=((1-x)*(roll[:-1]-value[:-1])[:,:,None,None]+x*(roll[1:]-value[1:])[:,:,None,None])/mesh.dz
    return (1-x)*lo+x*hi,dx,dz


def velocity(geometry,regular):
    """Input components (vr/r, vtheta/r, vz); return physical v and curl(v)."""
    mesh=geometry.mesh;regular=np.asarray(regular,float)
    if regular.shape!=(mesh.nr+1,mesh.nz,3) or not np.isfinite(regular).all():
        raise ValueError('Finite regular velocity components on dual nodes required')
    ur,_,urz=scalar(geometry,regular[...,0]);ut,utx,utz=scalar(geometry,regular[...,1])
    uz,uzx,_=scalar(geometry,regular[...,2]);r=geometry.r
    return np.stack((r*ur,r*ut,uz),axis=-1),np.stack((-r*utz,r*(urz-uzx),2*ut+r*r*utx),axis=-1)
