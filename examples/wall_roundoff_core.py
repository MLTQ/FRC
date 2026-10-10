"""Strict wall induction with balanced dual residual refinement."""
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from plasma_bench.inertial_induction import InductionRate
from plasma_bench.wall_induction_core import WallInductionCore


class WallRoundoffCore(WallInductionCore):
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
        solution = linear(rhs); self.refinement_trace = []; self.roundoff_details = []; self.roundoff_arrays = {}; norm = np.linalg.norm
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
            op_residual=operator(solution)-rhs
            physical_residual=self.dual_lift.T@difference
            dual_pre=self.dual_lift@solution;dual_post=geo.hessian@rate
            dual_gap=dual_pre-dual_post
            carried=geo.faraday@self.mass_inverse(model.inertia@self.mass_inverse(geo.faraday.T@dual_gap))
            metric_gap=self.lift.T@dual_pre-self.dual_lift.T@rate
            predicted=metric_gap+self.dual_lift.T@carried
            gap=op_residual+physical_residual
            self.roundoff_details.append(dict(refinement=refinement,
                dual_action_relative_difference=relative(dual_gap,dual_post),
                opposing_residual_gap_over_rhs=relative(gap,rhs),
                predicted_gap_over_rhs=relative(predicted,rhs),
                unexplained_fraction_of_gap=relative(gap-predicted,gap)))
            for key,value in dict(solution=solution,rhs=rhs,operator_residual=op_residual,
                physical_dual_residual=physical_residual,predicted_gap=predicted,dual_action_gap=dual_gap).items():
                self.roundoff_arrays[f'refinement_{refinement}_{key}']=value.copy()
            self.refinement_trace.append(dict(components))
            residual = max(components.values())
            if np.isfinite(residual) and residual <= self.tolerance:
                return InductionRate(rate, current, electric, iterations, residual, refinement, components)
            if refinement < 3:
                physical = self.dual_lift.T@difference
                # Algebraically identical equations have distinct rounding paths.
                # Balance their corrections only when the operator bound fails.
                correction = (.5*(physical+rhs-operator(solution))
                              if components['operator'] > self.tolerance else physical)
                solution += linear(correction)
        raise RuntimeError(f'Wall induction residuals {components} exceed tolerance after 3 refinements; no input modified')
