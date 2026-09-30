"""Gauss collocation with decoupled volume/boundary summation by parts."""
import numpy as np


class GaussSBP:
    def __init__(self,degree):
        if not isinstance(degree,(int,np.integer)) or degree<1:raise ValueError('Positive integer polynomial degree required')
        self.degree=degree;self.nodes,self.weights=np.polynomial.legendre.leggauss(degree+1)
        x=self.nodes;n=len(x);difference=x[:,None]-x[None,:];np.fill_diagonal(difference,1.)
        barycentric=1/np.prod(difference,axis=1)
        derivative=barycentric[None,:]/barycentric[:,None]/difference;np.fill_diagonal(derivative,0.)
        np.fill_diagonal(derivative,-derivative.sum(axis=1));self.derivative=derivative
        traces=barycentric[None,:]/(np.array([-1.,1.])[:,None]-x[None,:]);traces/=traces.sum(axis=1)[:,None]
        self.trace=traces;b=np.diag([-1.,1.]);q=self.weights[:,None]*derivative
        self.extended=np.block([[q-.5*traces.T@b@traces,.5*traces.T@b],[-.5*b@traces,.5*b]])
        self.lift=np.concatenate((np.eye(n),traces.T),axis=1)
        boundary=np.zeros_like(self.extended);boundary[-2:,-2:]=b
        if np.linalg.norm(self.extended+self.extended.T-boundary)>1e-10 or np.linalg.norm(self.extended@np.ones(n+2))>1e-10:
            raise RuntimeError('Gauss summation-by-parts identities failed')
