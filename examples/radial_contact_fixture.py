"""Sample a fixed independent radial contact mode on separate meshes."""
import numpy as np
from .strong_wall_reference import initialize
from .finite_mass_manufactured import cell_average


def fixture(spec,mode):
    if spec.get('kind','contact')!='contact':raise ValueError('Stationary contact fixture required')
    _,d,u,system,bias=initialize(dict(spec,kind='contact'));u[...,1:4]=0.
    target=mode.fields(d.r,d.z)['physical'];du=np.zeros_like(u);du[...,0]=target[...,0]
    np.divide(u[...,:1]*target[...,1:3],d.r[...,None],out=du[...,1:3],where=d.r[...,None]>0)
    du[-1,-1,:,1]=0.;du[...,3]=u[...,0]*target[...,3]
    du[...,4:]=u[...,:1]*target[...,4:6]+u[...,4:]/u[...,:1]*du[...,:1]
    mesh=system.geometry.mesh;psi=mode.fields(*mesh.corners)['psi']
    psi[[0,-1]]=0.;psi[:,-1]=psi[:,0]
    g=cell_average(mesh,lambda r,z:mode.fields(r,z)['physical'][...,7]/r,order=16)
    g=system.geometry.radial_projection@g
    return d,u,system,bias,du,system.geometry.pack(psi,g)
