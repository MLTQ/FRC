"""Same inertial induction equation on the finite-cylinder admissible field space."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.linalg import cho_factor
from scipy.sparse import block_diag,csr_matrix,eye,kron
from .vacuum_induction import VacuumInduction,_UniformFactor


class FiniteInduction(VacuumInduction):
    def __init__(self,model,*,tolerance=1e-10,max_iterations=1000):
        if not np.isfinite(tolerance) or not 0<tolerance<1 or type(max_iterations) is not int or max_iterations<1:
            raise ValueError('Positive tolerance and iteration budget required')
        self.model=model;self.tolerance=tolerance;self.max_iterations=max_iterations
        g=model.geometry;m=model.mesh;nodes=(m.nr+1)*(m.nz+1)
        active=np.r_[np.arange(m.nz+1,nodes),np.arange(nodes,nodes+m.nr*m.nz)]
        basis=block_diag((eye(m.nr*(m.nz+1)),kron(csr_matrix(g.radial_null),csr_matrix(g.axial_null))),format='csr')
        full=csr_matrix((np.ones(len(active)),(active,np.arange(len(active)))),shape=(g.hessian.shape[0],len(active)))@basis
        metric=(full.T@g.hessian@full).tocsc()
        self.lift=full.multiply((1/np.sqrt(metric.diagonal()))[None,:]).tocsr()
        self.dual_lift=(g.hessian@self.lift).tocsr()
        weak=(g.faraday.T@self.dual_lift).toarray()
        current=np.vstack([solve(weak[g.offsets[k]:g.offsets[k+1]]) for k,solve in enumerate(g.solves)])
        density=float(np.sum(g.weight*model.density)/np.sum(g.weight))
        metric=(self.lift.T@self.dual_lift).toarray()+model.reduced_mass/(E**2*density)*(weak.T@current)
        self.factor=_UniformFactor(cho_factor(.5*(metric+metric.T),lower=True))

    def separate_kernel(self,electric):
        g=self.model.geometry;m=self.model.mesh;r,t,z=g.split_electric(electric)
        r=r.reshape(m.nr,m.nz+1);z=z.reshape(m.nr+1,m.nz)
        rk=np.broadcast_to(np.mean(r,axis=1)[:,None],r.shape)
        zk=np.broadcast_to(z[:1],z.shape)
        small=np.r_[(r-rk).ravel(),t,(z-zk).ravel()]
        kernel=np.r_[rk.ravel(),np.zeros_like(t),zk.ravel()]
        return small,kernel
