"""Atomic full-equation implicit midpoint with a physical block chord solve."""
import numpy as np
from .plasma_coordinates import PlasmaCoordinates
from .plasma_preconditioner import PlasmaPreconditioner
from .regular_current import RegularCurrent
from .plasma_newton_direction import newton_direction


class ImplicitPlasmaMidpoint:
    def __init__(self, system, state, field, *, absolute_tolerance=1e-12,
                 relative_tolerance=1e-10, max_iterations=60):
        if not np.isfinite([absolute_tolerance, relative_tolerance]).all() or min(absolute_tolerance, relative_tolerance) <= 0:
            raise ValueError('Positive finite nonlinear residual tolerances required')
        if type(max_iterations) is not int or max_iterations < 1:
            raise ValueError('Positive nonlinear iteration budget required')
        self.system = system
        self.coordinates = PlasmaCoordinates(system, state, field)
        self.preconditioner = PlasmaPreconditioner(self.coordinates)
        self.atol = absolute_tolerance; self.rtol = relative_tolerance; self.max_iterations = max_iterations

    def step(self, state, field, dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite implicit time step required')
        c = self.coordinates; system = self.system
        state = np.asarray(state, float); system.fluid.fields(state)
        if field.axial_bias != c.bias:
            raise ValueError('Same prescribed constant harmonic bias required')
        initial_magnetic = field.packed_increment(system.geometry)
        self.preconditioner.factor(.5*dt)
        evaluations = 0; max_induction = 0.; max_refinement = 0; history = []; newton_steps = 0

        def residual(change):
            nonlocal evaluations, max_induction, max_refinement
            du, db = c.unpack(change)
            rate = system.derivative(state+.5*du, c.field(initial_magnetic+.5*db))
            evaluations += 1; max_induction = max(max_induction, rate.residual)
            max_refinement = max(max_refinement, rate.diagnostics['induction_residuals'].get('outer_refinements', 0))
            rhs = dt*c.pack(rate.fluid, rate.magnetic)
            error = change-rhs
            limit = self.atol+self.rtol*max(np.max(abs(change)), np.max(abs(rhs)))
            return error, float(np.max(abs(error))), float(limit)

        change = np.zeros(c.size)
        error, error_norm, limit = residual(change)
        for iteration in range(self.max_iterations+1):
            history.append(dict(iteration=iteration, residual=error_norm, limit=limit))
            if error_norm <= limit:
                du, db = c.unpack(change)
                updated = state+du; evolved = c.field(initial_magnetic+db)
                system.fluid.fields(updated)
                RegularCurrent(system.mapping, updated[..., 0], system.ion_mass)
                # Check the actual rounded returned states against the original
                # full time equation, not just an internal increment variable.
                actual = c.pack(updated-state, evolved.packed_increment(system.geometry)-initial_magnetic)
                _, stored_residual, stored_limit = residual(actual)
                if stored_residual > stored_limit:
                    raise RuntimeError('Rounded accepted-state midpoint residual failed; no input modified')
                return updated, evolved, dict(iterations=iteration, rhs_evaluations=evaluations,
                    newton_corrections=newton_steps,
                    residual=stored_residual, residual_limit=stored_limit, history=history,
                    max_induction_residual=max_induction, max_induction_refinements=max_refinement)
            if iteration == self.max_iterations:
                break
            accepted = False
            for method in ('chord', 'newton'):
                if method == 'chord':
                    correction = self.preconditioner.solve(-error)
                else:
                    correction, linear = newton_direction(residual, change, error, self.preconditioner.solve)
                    newton_steps += 1; history[-1].update(linear)
                for backtrack in range(13):
                    trial = change+2.**(-backtrack)*correction
                    try:
                        candidate_error, candidate_norm, candidate_limit = residual(trial)
                    except ValueError:
                        continue
                    if candidate_norm <= candidate_limit or candidate_norm < error_norm:
                        history[-1].update(backtracks=backtrack, correction=method)
                        change, error, error_norm, limit = trial, candidate_error, candidate_norm, candidate_limit
                        accepted = True; break
                if accepted: break
            if not accepted:
                raise RuntimeError(f'Implicit chord line search failed at residual {error_norm}; no input modified')
        raise RuntimeError(f'Implicit midpoint did not converge: residual {error_norm}, limit {limit}; no input modified')
