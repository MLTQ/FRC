"""Compatible full finite-cylinder fields with side/end magnetic vacuum energy."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import null_space
from scipy.sparse import block_diag,coo_matrix,csr_matrix,diags,kron,eye
from scipy.sparse.linalg import splu
from .compatible_geometry import CompatibleGeometry,gram
from .compatible_polynomials import bases
from .finite_cylinder_vacuum import FiniteCylinderVacuum


class FiniteCylinderGeometry(CompatibleGeometry):
    def __init__(self, mesh, *, degree=5,quadrature_order=12,vacuum_cells=8,modes=8):
        self.mesh=mesh;self.degree=degree;self.order=quadrature_order
        if type(quadrature_order) is not int or quadrature_order<degree+1:raise ValueError('Resolved physical quadrature required')
        q,w=np.polynomial.legendre.leggauss(quadrature_order)
        self.radial_points=(mesh.r[:,None]+.5*mesh.dr*q).ravel();self.axial_points=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
        self.wr=(np.pi*mesh.dr*self.radial_points.reshape(mesh.nr,-1)*w).ravel()
        self.wz=np.broadcast_to(.5*mesh.dz*w,(mesh.nz,len(q))).ravel()
        r=self.radial_points
        self.Nx,self.dNx,self.Ex,self.Dx=bases(.5*mesh.rf**2,.5*r*r,degree)
        self.Nz,self.dNz,self.Ez,self.Dz=bases(mesh.zf,self.axial_points,degree)
        self.radial_e=self.Ex.multiply(r[:,None]).tocsr()
        self.theta_e=(self.Nx[:,1:].multiply(1/r[:,None])@diags(mesh.rf[1:])).tocsr()
        self.radial_b=(self.Nx.multiply(1/r[:,None])@diags(np.r_[0.,np.ones(mesh.nr)])).tocsr()
        self.pairs=((self.radial_e,self.Nz),(self.theta_e,self.Nz),(self.Nx,self.Ez))
        self.mass=tuple(kron(gram(a,self.wr),gram(b,self.wz),format='csc') for a,b in self.pairs)
        self.solves=tuple(splu(m).solve for m in self.mass);self.offsets=np.cumsum([0]+[m.shape[0] for m in self.mass])
        hp=kron(gram(self.radial_b,self.wr),gram(self.dNz,self.wz))+kron(gram(self.dNx,self.wr),gram(self.Nz,self.wz))
        hg=kron(gram(self.radial_e,self.wr),gram(self.Ez,self.wz));self.plasma_hessian=block_diag((hp,hg),format='csc')/MU0
        self.vacuum=FiniteCylinderVacuum(mesh,degree=degree,quadrature_order=quadrature_order,vacuum_cells=vacuum_cells,modes=modes)
        # The vacuum's boundary order is side, lower end, upper end.
        self.boundary_indices=np.r_[mesh.nr*(mesh.nz+1)+np.arange(mesh.nz+1),
            np.arange(1,mesh.nr)*(mesh.nz+1),np.arange(1,mesh.nr)*(mesh.nz+1)+mesh.nz]
        self.boundary_map=csr_matrix((np.ones(len(self.boundary_indices)),(np.arange(len(self.boundary_indices)),self.boundary_indices)),
                                     shape=(len(self.boundary_indices),self.plasma_hessian.shape[0]))
        self.vacuum_hessian=(self.boundary_map.T@csr_matrix(self.vacuum.hessian)@self.boundary_map).tocsc()
        self.hessian=(self.plasma_hessian+self.vacuum_hessian).tocsc()
        self.radial_trace=bases(.5*mesh.rf**2,np.array([.5*mesh.radius**2]),degree)[2].toarray()
        self.axial_trace=bases(mesh.zf,mesh.zf[[0,-1]],degree)[2].toarray()
        def projection(mass,trace):
            dual=splu(mass).solve(trace.T);return np.eye(mass.shape[0])-dual@np.linalg.solve(trace@dual,trace)
        self.radial_projection=projection(gram(self.radial_e,self.wr),self.radial_trace)
        self.axial_projection=projection(gram(self.Ez,self.wz),self.axial_trace)
        self.radial_null=null_space(self.radial_trace);self.axial_null=null_space(self.axial_trace)
        nr,nz=mesh.nr,mesh.nz;nodes=(nr+1)*(nz+1)
        er_count=nr*(nz+1);t=kron(diags(-mesh.rf[1:]),eye(nz+1),format='csr')
        p=csr_matrix((nodes,self.offsets[-1]));rows=[];cols=[];values=[]
        for i in range(1,nr+1):
            for j in range(nz+1):rows.append(i*(nz+1)+j);cols.append(er_count+(i-1)*(nz+1)+j);values.append(-mesh.rf[i])
        for i in range(nr):
            for j in range(nz):
                for col,value in ((i*(nz+1)+j,1/mesh.dz),(i*(nz+1)+j+1,-1/mesh.dz),
                    (self.offsets[2]+i*nz+j,-1/mesh.radial_measure[i]),(self.offsets[2]+(i+1)*nz+j,1/mesh.radial_measure[i])):
                    rows.append(nodes+i*nz+j);cols.append(col);values.append(value)
        raw=coo_matrix((values,(rows,cols)),shape=(nodes+nr*nz,self.offsets[-1])).tocsr()
        self.projection=block_diag((eye(nodes),kron(csr_matrix(self.radial_projection),csr_matrix(self.axial_projection))),format='csr')
        self.faraday=(self.projection@raw).tocsr()
        self.r=np.broadcast_to(r.reshape(nr,1,len(q),1),(nr,nz,len(q),len(q)))
        self.z=np.broadcast_to(self.axial_points.reshape(1,nz,1,len(q)),self.r.shape)
        self.weight=self._layout(self.wr[:,None]*self.wz[None,:])

    def pack(self,psi,g):
        m=self.mesh;psi=np.asarray(psi);g=np.asarray(g)
        if psi.shape!=(m.nr+1,m.nz+1) or g.shape!=(m.nr,m.nz) or not np.isfinite(psi).all() or not np.isfinite(g).all():
            raise ValueError('Complete finite magnetic arrays required')
        if np.any(psi[0]!=0):raise ValueError('Regular zero axis flux required')
        for trace,value in ((self.radial_trace,g),(self.axial_trace,g.T)):
            relative=np.linalg.norm(trace@value)/max(np.linalg.norm(trace)*np.linalg.norm(value),np.finfo(float).tiny)
            if relative>256*np.finfo(float).eps:raise ValueError('Zero toroidal field at all vacuum interfaces required; no state repair')
        return np.r_[psi.ravel(),g.ravel()]

    def unpack(self,value):
        nodes=(self.mesh.nr+1)*(self.mesh.nz+1)
        return np.asarray(value[:nodes]).reshape(self.mesh.nr+1,self.mesh.nz+1),np.asarray(value[nodes:]).reshape(self.mesh.nr,self.mesh.nz)

    def magnetic_samples(self,psi,g):
        self.pack(psi,g)
        return np.stack((-self.tensor(self.radial_b,psi,self.dNz),self.tensor(self.radial_e,g,self.Ez),self.tensor(self.dNx,psi,self.Nz)),axis=-1)

    def gradient(self,potential):
        p=np.asarray(potential)
        if p.shape!=(self.mesh.nr+1,self.mesh.nz+1) or not np.isfinite(p).all():raise ValueError('Finite complete nodal potential required')
        return np.r_[(self.Dx@p).ravel(),np.zeros(self.mesh.nr*(self.mesh.nz+1)),np.asarray(p@self.Dz.T).ravel()]

    def project_curl(self,value):
        p,g=self.unpack(value)
        return np.r_[p.ravel(),(self.radial_projection@g@self.axial_projection.T).ravel()]

    def difference_curl(self,electric):
        m=self.mesh;r,t,z=self.split_electric(electric)
        r=r.reshape(m.nr,m.nz+1);t=t.reshape(m.nr,m.nz+1);z=z.reshape(m.nr+1,m.nz)
        psi=np.vstack((np.zeros((1,m.nz+1)),-m.rf[1:,None]*t))
        g=-np.diff(r,axis=1)/m.dz+np.diff(z,axis=0)/m.radial_measure[:,None]
        return self.project_curl(np.r_[psi.ravel(),g.ravel()])
