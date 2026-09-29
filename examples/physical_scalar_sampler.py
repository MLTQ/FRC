"""Independent physical sampling from root-product fields and primitive density."""
import numpy as np
from .compatible_field_norms import direct_basis
from .variable_current_norms import primitive_map


class PhysicalScalarSampler:
    def __init__(self,mesh,degree,order=12):
        self.mesh=mesh;self.order=order
        q,w=np.polynomial.legendre.leggauss(order)
        r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
        self.r=r[:,None];self.z=z[None]
        self.Nx,self.dNx,self.Ex=direct_basis(.5*mesh.rf**2,.5*r*r,degree)
        self.Nz,self.dNz,self.Ez=direct_basis(mesh.zf,z,degree,True)
        self.Rx=primitive_map(.5*mesh.rf**2,.5*r*r,degree)
        self.Rz=primitive_map(mesh.zf,z,degree,True)
        self.weight=(np.pi*mesh.dr*r.reshape(mesh.nr,-1)*w).ravel()[:,None]*np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()[None]

    def scalar(self,values):return self.Nx@values@self.Nz.T

    def density(self,values):return self.Rx@values@self.Rz.T

    def electric(self,value):
        mesh=self.mesh;cells=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
        er=value[:cells].reshape(mesh.nr,mesh.nz)
        et=np.pad(value[cells:cells+theta].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
        ez=value[cells+theta:].reshape(mesh.nr+1,mesh.nz)
        return np.stack((self.r*(self.Ex@er@self.Nz.T),
            (self.Nx@(mesh.rf[:,None]*et)@self.Nz.T)/self.r,self.Nx@ez@self.Ez.T),axis=-1)

    def magnetic(self,value):
        mesh=self.mesh;nodes=(mesh.nr+1)*mesh.nz
        p=value[:nodes].reshape(mesh.nr+1,mesh.nz);g=value[nodes:].reshape(mesh.nr,mesh.nz)
        return np.stack((-(self.Nx@p@self.dNz.T)/self.r,self.r*(self.Ex@g@self.Ez.T),
                         self.dNx@p@self.Nz.T),axis=-1)

    def cells(self,value):
        mesh=self.mesh
        return np.sum((self.weight*value).reshape(mesh.nr,self.order,mesh.nz,self.order),axis=(1,3))/mesh.volumes

    def metric(self,value,target):
        error=(value-target)**2;norm=target**2
        if error.ndim==3:error=error.sum(axis=-1);norm=norm.sum(axis=-1)
        total=float(np.sum(self.weight*error));denominator=float(np.sum(self.weight*norm))
        axis=self.order;axis_error=float(np.sum(self.weight[:axis]*error[:axis]))
        axis_norm=float(np.sum(self.weight[:axis]*norm[:axis]))
        return dict(relative_L2=float(np.sqrt(total/denominator)) if denominator>0 else None,
            rms_error=float(np.sqrt(total/self.mesh.volumes.sum())),sampled_maximum=float(np.sqrt(error.max())),
            axis_relative_L2=float(np.sqrt(axis_error/axis_norm)) if axis_norm>0 else None,
            axis_rms_error=float(np.sqrt(axis_error/self.weight[:axis].sum())))
