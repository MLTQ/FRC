"""Full-trace inertial induction including the eliminated vacuum field energy."""
from copy import copy
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.linalg import null_space, cho_factor, cho_solve
from scipy.sparse import block_diag, csr_matrix, eye, kron
from .split_inverse_wall_core import SplitInverseWallCore
from .inertial_induction import InductionRate
from .kernel_induction import separate


class VacuumInduction(SplitInverseWallCore):
    def __init__(self, model, *, tolerance=1e-10, max_iterations=1000):
        if not np.isfinite(tolerance) or not 0<tolerance<1 or type(max_iterations) is not int or max_iterations<1:
            raise ValueError('Positive finite tolerance and iteration budget required')
        self.model=model; self.tolerance=tolerance; self.max_iterations=max_iterations
        geo=model.geometry; m=model.mesh; nodes=(m.nr+1)*m.nz
        active=np.r_[np.arange(m.nz,nodes),np.arange(nodes,nodes+m.nr*m.nz)]
        radial=null_space(geo.wall_trace[None,:])
        basis=block_diag((eye(m.nr*m.nz),kron(csr_matrix(radial),eye(m.nz))),format='csr')
        full=csr_matrix((np.ones(len(active)),(active,np.arange(len(active)))),
                         shape=(geo.hessian.shape[0],len(active)))@basis
        metric=(full.T@geo.hessian@full).tocsc()
        self.lift=full.multiply((1/np.sqrt(metric.diagonal()))[None,:]).tocsr()
        self.dual_lift=(geo.hessian@self.lift).tocsr()
        weak=(geo.faraday.T@self.dual_lift).toarray()
        current=np.vstack([solve(weak[geo.offsets[k]:geo.offsets[k+1]]) for k,solve in enumerate(geo.solves)])
        density=float(np.sum(geo.weight*model.density)/np.sum(geo.weight))
        metric=(self.lift.T@self.dual_lift).toarray()+model.reduced_mass/(E**2*density)*(weak.T@current)
        self.uniform_factor=cho_factor(.5*(metric+metric.T),lower=True)
        # SplitInverseWallCore expects a .solve preconditioner. Only this inverse
        # uses uniform density; every physical inertia action uses actual n.
        self.factor=_UniformFactor(self.uniform_factor)

    def separate_kernel(self, electric):
        return separate(self.model.geometry,electric)

    def solve(self, forcing):
        model=self.model; geo=model.geometry; f=np.asarray(forcing,float)
        if f.shape != (model.mass.shape[0],) or not np.isfinite(f).all():
            raise ValueError('Finite full weak electric forcing required')
        electric0=self.mass_inverse(f); small,kernel=self.separate_kernel(electric0)
        weak=model.mass@small; relative=lambda a,b:float(np.linalg.norm(a)/max(np.linalg.norm(b),np.finfo(float).tiny))
        attempts=[]
        for refinement in range(5):
            inner=copy(self); inner.tolerance=self.tolerance*10.**(-refinement)
            try: result=SplitInverseWallCore.solve(inner,weak)
            except RuntimeError as exc:
                attempts.append(dict(tolerance=inner.tolerance,error=str(exc)));continue
            electric=result.electric+kernel; inertia=model.inertia@result.current
            curl=geo.difference_curl(electric)
            components=dict(result.residual_components,
                full_ohm_backward=relative(model.mass@electric-inertia-f,abs(model.mass@electric)+abs(inertia)+abs(f)),
                full_faraday_backward=relative(curl-result.magnetic,abs(geo.faraday)@abs(electric)+abs(result.magnetic)),
                kernel_reconstruction_backward=relative(electric0-small-kernel,abs(electric0)+abs(small)+abs(kernel)))
            residual=max(components.values()); attempts.append(dict(tolerance=inner.tolerance,residual=residual))
            if np.isfinite(residual) and residual<=self.tolerance:
                components.update(outer_refinements=refinement,inner_tolerance=inner.tolerance,
                    full_faraday_rate_relative=relative(curl-result.magnetic,result.magnetic))
                return InductionRate(result.magnetic,result.current,electric,result.iterations,residual,result.refinements,components)
        raise RuntimeError(f'Vacuum induction residual failure: {attempts}; no input modified')


class _UniformFactor:
    def __init__(self, factor):self.factor=factor
    def solve(self, value):return cho_solve(self.factor,value)
