"""Collocated endpoint/volume Gauss-Lobatto summation by parts."""
import numpy as np


class LobattoSBP:
    def __init__(self,degree):
        if not isinstance(degree,(int,np.integer)) or degree<2:raise ValueError('Integer Lobatto degree at least two required')
        p=np.polynomial.legendre.Legendre.basis(degree)
        roots=p.deriv().roots()
        if np.max(abs(np.imag(roots)))>1e-13:raise RuntimeError('Nonreal Lobatto nodes')
        self.nodes=np.r_[-1.,np.sort(np.real(roots)),1.];self.weights=2/(degree*(degree+1)*p(self.nodes)**2)
        difference=self.nodes[:,None]-self.nodes[None,:];np.fill_diagonal(difference,1.)
        barycentric=1/np.prod(difference,axis=1);d=barycentric[None,:]/barycentric[:,None]/difference;np.fill_diagonal(d,0.);np.fill_diagonal(d,-np.sum(d,axis=1))
        self.derivative=d;self.q=self.weights[:,None]*d
        boundary=np.diag(np.r_[-1.,np.zeros(degree-1),1.])
        if np.linalg.norm(self.q+self.q.T-boundary)>1e-10:raise RuntimeError('Lobatto SBP identity failed')
        # The regular radial/angular velocity at x=0 is the trace of the
        # degree-p-1 polynomial through the positive nodes, not another DOF.
        x=.5*(self.nodes[1:]+1);trace=np.ones(degree)
        for i in range(degree):
            for j in range(degree):
                if i!=j:trace[i]*=-x[j]/(x[i]-x[j])
        self.axis_trace=trace
