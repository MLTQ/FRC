"""Finite cylindrical material connectivity, independent of fluid dynamics."""
import numpy as np


class MaterialCylinder:
    def __init__(self,mesh,*,quadrature=6):
        self.nr=mesh.nr;self.nt=mesh.nz;self.radius=mesh.radius;self.quadrature=quadrature
        self.radial_grading=1.;r,z=mesh.corners
        self.fluid_reference=np.stack((r,z),axis=-1).reshape(-1,2)
        self.fluid_node_count=len(self.fluid_reference);self.fluid_cell_count=mesh.nr*mesh.nz
        self.ids=np.arange(self.fluid_node_count).reshape(mesh.nr+1,mesh.nz+1)
        self.cells=np.array([[self.ids[i,j],self.ids[i+1,j],self.ids[i+1,j+1],self.ids[i,j+1]] for i in range(mesh.nr) for j in range(mesh.nz)])
        self.axis=self.fluid_reference[:,0]==0
        self.boundary=np.r_[self.ids[-1],self.ids[:-1,0],self.ids[:-1,-1]]
        x,w=np.polynomial.legendre.leggauss(quadrature)
        xi,tau=np.meshgrid(.5*(x+1),.5*(x+1),indexing='ij');xi=xi.ravel();tau=tau.ravel()
        self.qweight=(.25*w[:,None]*w[None,:]).ravel()
        self.shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
        self.derivative=np.stack((np.stack((-(1-tau),1-tau,tau,-tau),axis=-1),np.stack((-(1-xi),-xi,xi,1-xi),axis=-1)),axis=-1)
