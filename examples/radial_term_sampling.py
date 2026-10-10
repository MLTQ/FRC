"""Bounded-memory independent sampling of archived physical term rates."""
from copy import copy
import numpy as np
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.mhd_mesh import MHDMesh


def radial_slice(sample,start,stop):
    result=copy(sample);rows=slice(start,stop)
    result.r=sample.r[rows];result.R=result.r[:,None];result.weight=sample.weight[rows]
    for key in ('nx','dx','ex','dex'):setattr(result,key,getattr(sample,key)[rows])
    result.radial={k:(v[0][rows],v[1][rows]) for k,v in sample.radial.items()}
    return result


def velocity_rate(sample,state,rate):
    n=sample.scalar(state[...,0]);nt=sample.scalar(rate[...,0]);vz=sample.scalar(state[...,3])/n
    vr=sample.R*sample.scalar((rate[...,1]-state[...,1]*rate[...,0]/state[...,0])/state[...,0],True)
    vt=sample.R*sample.scalar((rate[...,2]-state[...,2]*rate[...,0]/state[...,0])/state[...,0],True)
    return np.stack((vr,vt,(sample.scalar(rate[...,3])-vz*nt)/n),axis=-1)


def negative_curl(sample,electric):
    _,g=sample.electric(electric,True)
    return np.stack((g[...,1,2]-g[...,2,1],g[...,2,0]-g[...,0,2],g[...,0,1]-g[...,1,0]),axis=-1)


def raw_faraday(mesh,electric):
    n=mesh.nr*mesh.nz;t=(mesh.nr-1)*mesh.nz
    er=electric[:n].reshape(mesh.nr,mesh.nz);et=electric[n:n+t].reshape(mesh.nr-1,mesh.nz)
    ez=electric[n+t:].reshape(mesh.nr+1,mesh.nz);psi=np.zeros((mesh.nr+1,mesh.nz))
    psi[1:-1]=-mesh.rf[1:-1,None]*et
    g=(er-np.roll(er,-1,axis=1))/mesh.dz+np.diff(ez,axis=0)/mesh.radial_measure[:,None]
    return np.r_[psi.ravel(),g.ravel()]


def wall_projection(mesh,order=12):
    # The radial mass/trace projection is independent of the axial tensor factor.
    small=MHDMesh(mesh.radius,mesh.half_length,mesh.nr,8)
    return StrongWallGeometry(CompatibleGeometry(small,degree=5,quadrature_order=order)).radial_projection


def projected_faraday(mesh,electric,projection):
    result=raw_faraday(mesh,electric);count=(mesh.nr+1)*mesh.nz
    result[count:]=(projection@result[count:].reshape(mesh.nr,mesh.nz)).ravel()
    return result
