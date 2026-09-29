"""Higher-order cylindrical de Rham pairing with physical tensor-product masses."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import block_diag,coo_matrix,diags,kron
from scipy.sparse.linalg import splu
from .compatible_polynomials import bases


def gram(a,w):return (a.T@a.multiply(w[:,None])).tocsc()


class CompatibleGeometry:
    """Fixed psi at axis/wall, natural Btheta(R)=0, periodic z; no vacuum."""
    wall_btheta=0.

    def __init__(self,mesh,*,degree=5,quadrature_order=10):
        if type(quadrature_order) is not int or quadrature_order<degree+1:
            raise ValueError('Physical quadrature must resolve the selected polynomial degree')
        self.mesh=mesh;self.degree=degree;self.order=quadrature_order
        q,w=np.polynomial.legendre.leggauss(quadrature_order)
        self.radial_points=(mesh.r[:,None]+.5*mesh.dr*q).ravel()
        self.axial_points=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
        self.wr=(np.pi*mesh.dr*self.radial_points.reshape(mesh.nr,-1)*w).ravel()
        self.wz=np.broadcast_to(.5*mesh.dz*w,(mesh.nz,len(q))).ravel()
        self.Nx,self.dNx,self.Ex,self.Dx=bases(.5*mesh.rf**2,.5*self.radial_points**2,degree)
        self.Nz,self.dNz,self.Ez,self.Dz=bases(mesh.zf,self.axial_points,degree,periodic=True)
        radius=self.radial_points
        self.radial_e=self.Ex.multiply(radius[:,None]).tocsr()
        self.theta_e=(self.Nx[:,1:-1].multiply(1/radius[:,None])@diags(mesh.rf[1:-1])).tocsr()
        self.radial_b=self.Nx.multiply(1/radius[:,None]).tocsr()
        # Axis psi is constrained to zero; an inadmissible N_axis/r is not an
        # energy basis and must not introduce a quadrature-dependent singularity.
        self.radial_b=self.radial_b@diags(np.r_[0.,np.ones(mesh.nr)])
        self.pairs=((self.radial_e,self.Nz),(self.theta_e,self.Nz),(self.Nx,self.Ez))
        self.mass=tuple(kron(gram(a,self.wr),gram(b,self.wz),format='csc') for a,b in self.pairs)
        self.solves=tuple(splu(m).solve for m in self.mass)
        self.offsets=np.cumsum([0]+[m.shape[0] for m in self.mass])
        hpsi=kron(gram(self.radial_b,self.wr),gram(self.dNz,self.wz))+kron(gram(self.dNx,self.wr),gram(self.Nz,self.wz))
        hg=kron(gram(self.radial_e,self.wr),gram(self.Ez,self.wz))
        self.hessian=block_diag((hpsi,hg),format='csc')/MU0
        self.faraday=self._faraday()
        self.r=np.broadcast_to(radius.reshape(mesh.nr,1,len(q),1),(mesh.nr,mesh.nz,len(q),len(q)))
        self.z=np.broadcast_to(self.axial_points.reshape(1,mesh.nz,1,len(q)),self.r.shape)
        self.weight=self._layout(self.wr[:,None]*self.wz[None])

    def _faraday(self):
        mesh=self.mesh;nr,nz=mesh.nr,mesh.nz;nodes=(nr+1)*nz
        rows=[];cols=[];values=[]
        for i in range(nr):
            for j in range(nz):
                for col,value in ((i*nz+j,1/mesh.dz),(i*nz+(j+1)%nz,-1/mesh.dz),
                    (self.offsets[2]+i*nz+j,-1/mesh.radial_measure[i]),(self.offsets[2]+(i+1)*nz+j,1/mesh.radial_measure[i])):
                    rows.append(nodes+i*nz+j);cols.append(col);values.append(value)
        for i in range(1,nr):
            for j in range(nz):rows.append(i*nz+j);cols.append(self.offsets[1]+(i-1)*nz+j);values.append(-mesh.rf[i])
        return coo_matrix((values,(rows,cols)),shape=(nodes+nr*nz,self.offsets[-1])).tocsr()

    def _layout(self,value):return np.asarray(value).reshape(self.mesh.nr,self.order,self.mesh.nz,self.order).transpose(0,2,1,3)

    def _unlayout(self,value):return np.asarray(value).transpose(0,2,1,3).reshape(self.mesh.nr*self.order,self.mesh.nz*self.order)

    def tensor(self,a,value,b):return self._layout((a@value)@b.T)

    def pack(self,psi,g):
        if np.shape(psi)!=(self.mesh.nr+1,self.mesh.nz+1) or np.shape(g)!=(self.mesh.nr,self.mesh.nz):raise ValueError('Complete periodic magnetic shapes required')
        if not np.isfinite(psi).all() or not np.isfinite(g).all() or np.any(psi[0]!=0) or np.any(psi[:,0]!=psi[:,-1]):raise ValueError('Finite regular-axis periodic flux required')
        return np.r_[psi[:,:-1].ravel(),g.ravel()]

    def unpack(self,value):
        nodes=(self.mesh.nr+1)*self.mesh.nz
        psi=np.asarray(value[:nodes]).reshape(self.mesh.nr+1,self.mesh.nz)
        return np.c_[psi,psi[:,0]],np.asarray(value[nodes:]).reshape(self.mesh.nr,self.mesh.nz)

    def split_electric(self,value):return tuple(value[self.offsets[k]:self.offsets[k+1]] for k in range(3))

    def sample_electric(self,value):
        result=[]
        for (a,b),part in zip(self.pairs,self.split_electric(value)):
            result.append(self.tensor(a,part.reshape(a.shape[1],b.shape[1]),b))
        return np.stack(result,axis=-1)

    def project_electric(self,values):
        if np.shape(values)!=self.r.shape+(3,) or not np.isfinite(values).all():raise ValueError('Finite physical vector at common quadrature required')
        result=[]
        for k,((a,b),solve) in enumerate(zip(self.pairs,self.solves)):
            weighted=self._unlayout(values[...,k])*self.wr[:,None]*self.wz[None]
            result.append(solve(np.asarray((a.T@weighted)@b).ravel()))
        return np.concatenate(result)

    def magnetic_samples(self,psi,g):
        self.pack(psi,g)
        p=psi[:,:-1]
        return np.stack((-self.tensor(self.radial_b,p,self.dNz),
                         self.tensor(self.radial_e,g,self.Ez),self.tensor(self.dNx,p,self.Nz)),axis=-1)

    def gradient(self,potential):
        p=np.asarray(potential)
        if p.shape!=(self.mesh.nr+1,self.mesh.nz) or not np.isfinite(p).all():raise ValueError('Finite full nodal potential required')
        return np.r_[(self.Dx@p).ravel(),np.zeros((self.mesh.nr-1)*self.mesh.nz),np.asarray(p@self.Dz.T).ravel()]
