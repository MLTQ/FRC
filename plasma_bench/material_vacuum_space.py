"""Full material field/current coordinates with a matched exterior flux trace."""
import numpy as np
from scipy.sparse import csr_matrix,hstack,kron
from .compatible_polynomials import bases
from .material_current_space import MaterialCurrentSpace


class MaterialVacuumSpace:
    def __init__(self,mesh,*,degree=3,quadrature=6):
        base=MaterialCurrentSpace(mesh,degree=degree,quadrature=quadrature)
        self.mesh=mesh;self.degree=degree;self.order=quadrature;self.reference_r=base.reference_r;self.reference_z=base.reference_z
        q=quadrature;r=base.reference_r.reshape(mesh.nr,mesh.nz,q,q)[:,0,:,0].ravel();z=base.reference_z.reshape(mesh.nr,mesh.nz,q,q)[0,:,0,:].ravel()
        nx,dx,ex,_=bases(.5*mesh.rf**2,.5*r*r,degree);nz,dz,ez,_=bases(mesh.zf,z,degree)
        order=np.arange(len(r)*len(z)).reshape(mesh.nr,q,mesh.nz,q).transpose(0,2,1,3).ravel()
        # psi is the same nodal polynomial as the exterior; psi(axis)=0.
        flux=kron(nx[:,1:],nz,format='csr')[order]
        etheta=flux.multiply((1/self.reference_r)[:,None]).tocsr()
        br=-kron(nx[:,1:],dz,format='csr')[order].multiply((1/self.reference_r)[:,None])
        bz=kron(dx[:,1:],nz,format='csr')[order]
        bt=kron(ex.multiply(r[:,None]),ez,format='csr')[order]
        ne=base.electric_count;nt=bt.shape[1];npoloidal=flux.shape[1];points=len(order)
        zero_p=csr_matrix((points,npoloidal));zero_e=csr_matrix((points,ne));zero_t=csr_matrix((points,nt))
        self.electric=(hstack((base.er,zero_p),format='csr'),hstack((zero_e,etheta),format='csr'),hstack((base.ez,zero_p),format='csr'))
        self.curl=(hstack((zero_e,br),format='csr'),hstack((base.curl,zero_p),format='csr'),hstack((zero_e,bz),format='csr'))
        self.magnetic=(hstack((zero_t,br),format='csr'),hstack((bt,zero_p),format='csr'),hstack((zero_t,bz),format='csr'))
        self.flux=hstack((zero_t,flux),format='csr');self.gradient=np.vstack((base.gradient,np.zeros((npoloidal,base.gradient.shape[1]))))
        self.toroidal_magnetic_count=nt;self.meridional_electric_count=ne
        self.electric_count=ne+npoloidal;self.magnetic_count=nt+npoloidal
        boundary_nodes=np.r_[np.arange(mesh.nr*(mesh.nz+1),(mesh.nr+1)*(mesh.nz+1)),np.arange(1,mesh.nr)*(mesh.nz+1),np.arange(1,mesh.nr)*(mesh.nz+1)+mesh.nz]
        self.vacuum_trace=csr_matrix((np.ones(len(boundary_nodes)),(np.arange(len(boundary_nodes)),nt+boundary_nodes-(mesh.nz+1))),shape=(len(boundary_nodes),self.magnetic_count))
