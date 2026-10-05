"""Full Gauss equations with accepted-residual secant proposals."""
import numpy as np
from .modal_wall_gauss import ModalWallGauss
from .regular_current import RegularCurrent
from .plasma_newton_direction import newton_direction
from .gauss_secant import GaussSecant


class AcceleratedModalGauss(ModalWallGauss):
    def step(self, state, field, dt):
        if not np.isfinite(dt) or dt <= 0: raise ValueError('Positive finite implicit time step required')
        c = self.coordinates; system = self.system; state = np.asarray(state, float)
        system.fluid.fields(state)
        if field.axial_bias != c.bias: raise ValueError('Same prescribed constant harmonic bias required')
        initial_magnetic = field.packed_increment(system.geometry)
        self.inverse.factor(dt*(.25+1j*np.sqrt(3.)/12))
        evaluations = 0; max_induction = 0.; max_refinement = 0; history = []; newton_steps = 0

        def residual(changes):
            nonlocal evaluations, max_induction, max_refinement
            rates = []; actual = []
            for k, change in enumerate(changes):
                du, db = c.unpack(change)
                stage = state+du; stage_magnetic = initial_magnetic+db
                actual.append(c.pack(stage-state, stage_magnetic-initial_magnetic))
                if k == 1 and np.array_equal(changes[0], changes[1]):
                    rates.append(rates[0]); continue
                rate = system.derivative(stage, c.field(stage_magnetic))
                evaluations += 1; max_induction = max(max_induction, rate.residual)
                max_refinement = max(max_refinement, rate.diagnostics['induction_residuals'].get('outer_refinements', 0))
                rates.append(c.pack(rate.fluid, rate.magnetic))
            rhs = dt*(self.a@np.asarray(rates))
            error = np.asarray(actual)-rhs
            limit = self.atol+self.rtol*max(np.max(abs(changes)), np.max(abs(rhs)))
            return error, float(np.max(abs(error))), float(limit), np.asarray(rates)

        secant = GaussSecant(self.inverse.stages)
        accelerated_steps = 0
        changes = np.zeros((2, c.size))
        error, error_norm, limit, rates = residual(changes)
        for iteration in range(self.max_iterations+1):
            history.append(dict(iteration=iteration, residual=error_norm, limit=limit))
            if error_norm <= limit:
                change = .5*dt*(rates[0]+rates[1]); du, db = c.unpack(change)
                updated = state+du; evolved = c.field(initial_magnetic+db)
                system.fluid.fields(updated); RegularCurrent(system.mapping, updated[..., 0], system.ion_mass)
                endpoint_error = float(np.max(abs(c.pack(updated-state,
                    evolved.packed_increment(system.geometry)-initial_magnetic)-change)))
                if endpoint_error > self.atol+self.rtol*np.max(abs(change)):
                    raise RuntimeError('Rounded accepted-state Gauss residual failed; no input modified')
                return updated, evolved, dict(iterations=iteration, rhs_evaluations=evaluations,
                    newton_corrections=newton_steps, accelerated_corrections=accelerated_steps,
                    residual=max(error_norm, endpoint_error), residual_limit=limit, history=history,
                    endpoint_residual=endpoint_error, stage_changes=changes.copy(),
                    max_induction_residual=max_induction, max_induction_refinements=max_refinement)
            if iteration == self.max_iterations: break
            accepted = False
            methods = ('secant', 'chord', 'newton') if secant.images else ('chord', 'newton')
            for method in methods:
                if method == 'secant':
                    correction, rank = secant.correction(error)
                    history[-1]['secant_rank'] = rank
                elif method == 'chord':
                    correction = self.inverse.stages(-error)
                else:
                    correction, linear = newton_direction(residual, changes, error, self.inverse.stages)
                    newton_steps += 1; history[-1].update(linear)
                for backtrack in range(13):
                    trial = changes+2.**(-backtrack)*correction
                    try:
                        candidate_error, candidate_norm, candidate_limit, candidate_rates = residual(trial)
                    except ValueError:
                        continue
                    if candidate_norm <= candidate_limit or candidate_norm < error_norm:
                        history[-1].update(backtracks=backtrack, correction=method)
                        secant.append(trial-changes, candidate_error-error)
                        accelerated_steps += int(method == 'secant')
                        changes, error, error_norm, limit, rates = trial, candidate_error, candidate_norm, candidate_limit, candidate_rates
                        accepted = True; break
                if accepted: break
            if not accepted:
                raise RuntimeError(f'Implicit Gauss line search failed at residual {error_norm}; no input modified')
        raise RuntimeError(f'Implicit Gauss did not converge: residual {error_norm}, limit {limit}; no input modified')
