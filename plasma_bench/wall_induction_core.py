"""Positive inertial induction with every Krylov vector in the wall space."""
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import block_diag, csr_matrix, eye, kron
from scipy.sparse.linalg import LinearOperator, cg, splu
from .inertial_induction import InductionRate


class WallInductionCore:
    def __init__(self, model, *, tolerance=2e-12, max_iterations=1000):
        if not np.isfinite(tolerance) or not 0 < tolerance < 1 or type(max_iterations) is not int or max_iterations < 1:
            raise ValueError('Positive finite tolerance below one and iteration budget required')
        self.model = model; self.tolerance = tolerance; self.max_iterations = max_iterations
        geo = model.geometry; mesh = model.mesh; nodes = (mesh.nr+1)*mesh.nz
        self.active = np.r_[np.arange(mesh.nz, mesh.nr*mesh.nz), np.arange(nodes, nodes+mesh.nr*mesh.nz)]
        radial = null_space(geo.wall_trace[None, :])
        basis = block_diag((eye((mesh.nr-1)*mesh.nz), kron(csr_matrix(radial), eye(mesh.nz))), format='csr')
        full = csr_matrix((np.ones(len(self.active)), (self.active, np.arange(len(self.active)))),
                          shape=(geo.hessian.shape[0], len(self.active)))@basis
        metric = (full.T@geo.hessian@full).tocsc()
        scale = 1/np.sqrt(metric.diagonal())
        self.lift = full.multiply(scale[None, :]).tocsr()
        self.dual_lift = (geo.hessian@self.lift).tocsr()
        self.factor = splu((self.lift.T@self.dual_lift).tocsc())

    def mass_inverse(self, value):
        geo = self.model.geometry
        return np.concatenate([solve(v) for solve, v in zip(geo.solves, geo.split_electric(value))])

    def solve(self, forcing):
        model = self.model; geo = model.geometry; f = np.asarray(forcing, float)
        if f.shape != (model.mass.shape[0],) or not np.isfinite(f).all():
            raise ValueError('Finite complete weak forcing required')
        def operator(value):
            dual = self.dual_lift@value
            current = self.mass_inverse(geo.faraday.T@dual)
            correction = geo.faraday@self.mass_inverse(model.inertia@current)
            return self.lift.T@dual+self.dual_lift.T@correction
        size = self.lift.shape[1]
        system = LinearOperator((size, size), matvec=operator, dtype=float)
        preconditioner = LinearOperator((size, size), matvec=self.factor.solve, dtype=float)
        rhs = self.dual_lift.T@(geo.faraday@self.mass_inverse(f)); iterations = 0
        def count(_):
            nonlocal iterations
            iterations += 1
        def linear(right):
            result, info = cg(system, right, M=preconditioner, rtol=.1*self.tolerance, atol=0.,
                              maxiter=self.max_iterations, callback=count)
            if info != 0: raise RuntimeError(f'Wall induction did not converge: {info}; no input modified')
            return result
        solution = linear(rhs); norm = np.linalg.norm
        relative = lambda a, b: float(norm(a)/max(norm(b), np.finfo(float).tiny))
        for refinement in range(4):
            rate = self.lift@solution
            current = -self.mass_inverse(geo.faraday.T@(geo.hessian@rate))
            electric = self.mass_inverse(f+model.inertia@current)
            ohm = model.mass@electric-model.inertia@current-f
            ampere = model.mass@current+geo.faraday.T@(geo.hessian@rate)
            difference = geo.faraday@electric-rate
            components = dict(operator=relative(operator(solution)-rhs, rhs),
                faraday_dual=relative(self.dual_lift.T@difference, rhs),
                ohm=relative(ohm, abs(f)+abs(model.inertia@current)),
                ampere=relative(ampere, abs(model.mass@current)+abs(geo.faraday.T@(geo.hessian@rate))))
            field_norm = np.sqrt(rate@(geo.hessian@rate))
            components['faraday_physical'] = float(np.sqrt(difference@(geo.hessian@difference))/max(field_norm, np.finfo(float).tiny))
            residual = max(components.values())
            if np.isfinite(residual) and residual <= self.tolerance:
                return InductionRate(rate, current, electric, iterations, residual, refinement, components)
            if refinement < 3: solution += linear(self.dual_lift.T@difference)
        raise RuntimeError(f'Wall induction residuals {components} exceed tolerance after 3 refinements; no input modified')
