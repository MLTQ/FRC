"""Physical number-weighted scalar metric on compatible continuous nodal fields."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import kron
from scipy.sparse.linalg import LinearOperator,cg,splu
from .compatible_geometry import gram


@dataclass(frozen=True)
class ScalarSolution:
    values: np.ndarray
    residual: float
    iterations: int


class EntropyScalar:
    """Uses the caller's immutable conservative density; no second density state."""
    def __init__(self,density_map,*,tolerance=2e-12,max_iterations=500):
        if not np.isfinite(tolerance) or not 0<tolerance<1 or type(max_iterations) is not int or max_iterations<1:
            raise ValueError('Positive finite scalar tolerance and iteration budget required')
        self.density_map=density_map;self.geometry=density_map.geometry;self.mesh=density_map.mesh
        self.tolerance=tolerance;self.max_iterations=max_iterations
        geo=self.geometry;self.shape=(self.mesh.nr+1,self.mesh.nz);self.size=np.prod(self.shape)
        self.n=density_map.samples
        self.scale=float(np.sum(geo.weight*self.n)/self.mesh.volumes.sum())
        self.normalized_density=self.n/self.scale
        mass=kron(gram(geo.Nx,geo.wr),gram(geo.Nz,geo.wz),format='csc')
        self.factor=splu(mass)
        self.operator=LinearOperator((self.size,self.size),matvec=self._operator,dtype=float)
        self.preconditioner=LinearOperator(self.operator.shape,matvec=self.factor.solve,dtype=float)

    def validate(self,values):
        values=np.asarray(values,float)
        if values.shape!=self.shape or not np.isfinite(values).all():
            raise ValueError('Finite periodic nodal scalar including axis and wall required')
        return values

    def sample(self,values):
        geo=self.geometry
        return geo.tensor(geo.Nx,self.validate(values),geo.Nz)

    def gradient(self,values):
        # Incidence takes differences before interpolation: constant entropy has
        # exactly zero gradient, including at the regular cylindrical axis.
        return self.geometry.sample_electric(self.geometry.gradient(self.validate(values)))

    def adjoint(self,values):
        geo=self.geometry;values=np.asarray(values,float)
        if values.shape!=geo.r.shape or not np.isfinite(values).all():
            raise ValueError('Finite scalar at common physical quadrature required')
        return np.asarray((geo.Nx.T@geo._unlayout(geo.weight*values))@geo.Nz)

    def _operator(self,vector):
        return self.adjoint(self.normalized_density*self.sample(np.asarray(vector).reshape(self.shape))).ravel()

    def mass_action(self,values):
        return self.scale*self._operator(self.validate(values).ravel()).reshape(self.shape)

    def solve(self,physical_load):
        rhs=self.validate(physical_load).ravel()/self.scale;iterations=0
        if not np.any(rhs):return ScalarSolution(np.zeros(self.shape),0.,0)
        def count(_):
            nonlocal iterations
            iterations+=1
        result,info=cg(self.operator,rhs,M=self.preconditioner,rtol=.1*self.tolerance,
                       atol=0.,maxiter=self.max_iterations,callback=count)
        residual=float(np.linalg.norm(self.operator@result-rhs)/np.linalg.norm(rhs))
        if info or not np.isfinite(result).all() or residual>self.tolerance:
            raise RuntimeError(f'Physical scalar solve failed: info={info}, residual={residual}; no state modified')
        return ScalarSolution(result.reshape(self.shape),residual,iterations)

    def project(self,values):
        """Number-weighted physical L2 projection, preserving an explicit constant."""
        values=np.asarray(values,float)
        if values.shape!=self.n.shape or not np.isfinite(values).all():
            raise ValueError('Finite scalar at common physical quadrature required')
        anchor=float(values.flat[0])
        result=self.solve(self.adjoint(self.n*(values-anchor)))
        return ScalarSolution(result.values+anchor,result.residual,result.iterations)
