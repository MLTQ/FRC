"""Finite compatible toroidal magnetic sector for a moving material cylinder."""
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import csr_matrix,diags,kron,hstack,eye
from .compatible_polynomials import bases


class MaterialCurrentSpace:
    def __init__(self,mesh,*,degree=3,quadrature=6):
        self.mesh=mesh;self.degree=degree;self.order=quadrature
        q,w=np.polynomial.legendre.leggauss(quadrature)
        r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
        nx,dx,ex,Dx=bases(.5*mesh.rf**2,.5*r*r,degree)
        nz,dz,ez,Dz=bases(mesh.zf,z,degree)
        order=np.arange(len(r)*len(z)).reshape(mesh.nr,quadrature,mesh.nz,quadrature).transpose(0,2,1,3).ravel()
        radial=kron(ex.multiply(r[:,None]),nz,format='csr')[order]
        axial=kron(nx,ez,format='csr')[order]
        self.er=hstack((radial,csr_matrix((len(order),axial.shape[1]))),format='csr')
        self.ez=hstack((csr_matrix((len(order),radial.shape[1])),axial),format='csr')
        raw_b=kron(ex.multiply(r[:,None]),ez,format='csr')[order]
        self.incidence=hstack((-kron(eye(mesh.nr),Dz),kron(Dx,eye(mesh.nz))),format='csr')
        self.curl=-raw_b@self.incidence
        rt=bases(.5*mesh.rf**2,np.array([.5*mesh.radius**2]),degree)[2].toarray()
        zt=bases(mesh.zf,mesh.zf[[0,-1]],degree)[2].toarray()
        self.null=kron(csr_matrix(null_space(rt)),csr_matrix(null_space(zt)),format='csr')
        self.magnetic=(raw_b@self.null).tocsr();self.electric_count=self.er.shape[1];self.magnetic_count=self.null.shape[1]
        self.gradient=np.vstack(((kron(Dx,eye(mesh.nz+1))).toarray(),(kron(eye(mesh.nr+1),Dz)).toarray()))
        self.reference_r=np.broadcast_to(r.reshape(mesh.nr,1,quadrature,1),(mesh.nr,mesh.nz,quadrature,quadrature)).ravel()
        self.reference_z=np.broadcast_to(z.reshape(1,mesh.nz,1,quadrature),(mesh.nr,mesh.nz,quadrature,quadrature)).ravel()
