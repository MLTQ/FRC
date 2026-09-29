"""General positive magnetic-space solve for instantaneous inertial induction."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator, cg, splu


@dataclass(frozen=True)
class InductionRate:
    magnetic: np.ndarray
    current: np.ndarray
    electric: np.ndarray
    iterations: int
    residual: float
    refinements: int
    residual_components: dict


class InertialInduction:
    """Solve M E=I J_t+f, M J_t=-C.T H B_t, B_t=C E at fixed instantaneous n."""
    def __init__(self,model,*,tolerance=2e-12,max_iterations=1000):
        if not np.isfinite(tolerance) or not 0<tolerance<1 or type(max_iterations) is not int or max_iterations<1:
            raise ValueError('Positive finite tolerance below one and iteration budget required')
        self.model=model;self.tolerance=tolerance;self.max_iterations=max_iterations
        mesh=model.mesh;nodes=(mesh.nr+1)*mesh.nz
        self.active=np.r_[np.arange(mesh.nz,mesh.nr*mesh.nz),np.arange(nodes,nodes+mesh.nr*mesh.nz)]
        self.curl=model.geometry.faraday[self.active]
        self.hessian=model.hessian[self.active][:,self.active].tocsc()
        self.scale=1/np.sqrt(self.hessian.diagonal())
        diagonal=diags(self.scale)
        self.scaled_hessian=(diagonal@self.hessian@diagonal).tocsc()
        self.factor=splu(self.scaled_hessian)

    def mass_inverse(self,v):
        geo=self.model.geometry
        return np.concatenate([solve(a) for solve,a in zip(geo.solves,geo.split_electric(v))])

    def solve(self,forcing):
        model=self.model;f=np.asarray(forcing,float)
        if f.shape!=(model.mass.shape[0],) or not np.isfinite(f).all():
            raise ValueError('Finite weak electric forcing with the complete electric shape required')
        def operator(y):
            magnetic=self.scale*y;dual=self.hessian@magnetic
            current=self.mass_inverse(self.curl.T@dual)
            correction=self.curl@self.mass_inverse(model.inertia@current)
            return self.scale*(dual+self.hessian@correction)
        size=len(self.active);system=LinearOperator((size,size),matvec=operator,dtype=float)
        preconditioner=LinearOperator((size,size),matvec=self.factor.solve,dtype=float)
        rhs=self.scale*(self.hessian@(self.curl@self.mass_inverse(f)));iterations=0
        def count(_):
            nonlocal iterations
            iterations+=1
        def linear(right):
            result,info=cg(system,right,M=preconditioner,rtol=.1*self.tolerance,atol=0.,
                           maxiter=self.max_iterations,callback=count)
            if info!=0:raise RuntimeError(f'Inertial induction did not converge: {info}; no input modified')
            return result
        solution=linear(rhs)
        norm=lambda v:np.linalg.norm(v)
        relative=lambda a,b:float(norm(a)/max(norm(b),np.finfo(float).tiny))
        for refinement in range(4):
            rate=np.zeros(model.geometry.faraday.shape[0]);rate[self.active]=self.scale*solution
            current=-self.mass_inverse(self.curl.T@(self.hessian@rate[self.active]))
            electric=self.mass_inverse(f+model.inertia@current)
            ohm=model.mass@electric-model.inertia@current-f
            ampere=model.mass@current+model.geometry.faraday.T@(model.hessian@rate)
            difference=(model.geometry.faraday@electric-rate)[self.active]
            components=dict(operator=relative(operator(solution)-rhs,rhs),
                faraday_dual=relative(self.scale*(self.hessian@difference),rhs),
                ohm=relative(ohm,abs(f)+abs(model.inertia@current)),
                ampere=relative(ampere,abs(model.mass@current)+abs(model.geometry.faraday.T@(model.hessian@rate))))
            field_norm=np.sqrt(rate[self.active]@(self.hessian@rate[self.active]))
            components['faraday_physical']=float(np.sqrt(difference@(self.hessian@difference))/max(field_norm,np.finfo(float).tiny))
            residual=max(components.values())
            if np.isfinite(residual) and residual<=self.tolerance:
                return InductionRate(rate,current,electric,iterations,residual,refinement,components)
            if refinement<3:
                solution+=linear(self.scale*(self.hessian@difference))
        raise RuntimeError(f'Inertial induction residuals {components} exceed tolerance after 3 refinements; no input modified')
