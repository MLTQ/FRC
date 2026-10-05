"""Full nonlinear collocation step with a declared final stage and weights."""
import numpy as np
from .plasma_newton_direction import newton_direction
from .regular_current import RegularCurrent


def step(solver, state, field, dt):
    if not np.isfinite(dt) or dt <= 0: raise ValueError('Positive finite implicit time step required')
    c = solver.coordinates; system = solver.system; state = np.asarray(state, float)
    system.fluid.fields(state)
    if field.axial_bias != c.bias: raise ValueError('Same prescribed constant harmonic bias required')
    initial_magnetic = field.packed_increment(system.geometry)
    solver.inverse.factor(dt*solver.shift)
    evaluations = 0; max_induction = 0.; max_refinement = 0; history = []; newton_steps = 0

    def residual(changes):
        nonlocal evaluations, max_induction, max_refinement
        rates = []; actual = []
        for k, change in enumerate(changes):
            du, db = c.unpack(change)
            stage = state+du; magnetic = initial_magnetic+db
            actual.append(c.pack(stage-state, magnetic-initial_magnetic))
            if k == 1 and np.array_equal(changes[0], changes[1]):
                rates.append(rates[0]); continue
            rate = system.derivative(stage, c.field(magnetic))
            evaluations += 1; max_induction = max(max_induction, rate.residual)
            max_refinement = max(max_refinement, rate.diagnostics['induction_residuals'].get('outer_refinements', 0))
            rates.append(c.pack(rate.fluid, rate.magnetic))
        rhs = dt*(solver.a@np.asarray(rates))
        error = np.asarray(actual)-rhs
        limit = solver.atol+solver.rtol*max(np.max(abs(changes)), np.max(abs(rhs)))
        return error, float(np.max(abs(error))), float(limit), np.asarray(rates)

    changes = np.zeros((2, c.size)); error, error_norm, limit, rates = residual(changes)
    for iteration in range(solver.max_iterations+1):
        history.append(dict(iteration=iteration, residual=error_norm, limit=limit))
        if error_norm <= limit:
            du, db = c.unpack(changes[solver.endpoint_stage])
            updated = state+du; evolved = c.field(initial_magnetic+db)
            system.fluid.fields(updated); RegularCurrent(system.mapping, updated[..., 0], system.ion_mass)
            rhs = dt*(solver.weights@rates)
            endpoint_error = float(np.max(abs(c.pack(updated-state,
                evolved.packed_increment(system.geometry)-initial_magnetic)-rhs)))
            endpoint_limit = solver.atol+solver.rtol*max(np.max(abs(changes[solver.endpoint_stage])), np.max(abs(rhs)))
            if endpoint_error > endpoint_limit:
                raise RuntimeError('Returned collocation endpoint equation failed; no input modified')
            return updated, evolved, dict(iterations=iteration, rhs_evaluations=evaluations,
                newton_corrections=newton_steps, residual=max(error_norm, endpoint_error),
                residual_limit=limit, history=history, endpoint_residual=endpoint_error,
                stage_changes=changes.copy(), max_induction_residual=max_induction,
                max_induction_refinements=max_refinement)
        if iteration == solver.max_iterations: break
        accepted = False
        for method in ('chord', 'newton'):
            if method == 'chord': correction = solver.inverse.stages(-error)
            else:
                correction, linear = newton_direction(residual, changes, error, solver.inverse.stages)
                newton_steps += 1; history[-1].update(linear)
            for backtrack in range(13):
                trial = changes+2.**(-backtrack)*correction
                try: candidate_error, candidate_norm, candidate_limit, candidate_rates = residual(trial)
                except ValueError: continue
                if candidate_norm <= candidate_limit or candidate_norm < error_norm:
                    history[-1].update(backtracks=backtrack, correction=method)
                    changes, error, error_norm, limit, rates = trial, candidate_error, candidate_norm, candidate_limit, candidate_rates
                    accepted = True; break
            if accepted: break
        if not accepted: raise RuntimeError(f'Collocation line search failed at residual {error_norm}; no input modified')
    raise RuntimeError(f'Collocation did not converge: residual {error_norm}, limit {limit}; no input modified')
