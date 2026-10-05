"""Complete unforced implicit plasma trajectories and replayable raw states."""
import time
import numpy as np
from plasma_bench.strong_wall_gauss import StrongWallGauss
from plasma_bench.regular_current import RegularCurrent
from .strong_wall_reference import initialize
from .regular_plasma_evolution import field_norms, grid_change


def run_case(spec, progress=None):
    start = time.monotonic()
    ref, d, initial, system, initial_field = initialize(spec)
    method = spec.get('method', 'gauss')
    classes = {'gauss': StrongWallGauss}
    if method not in classes: raise ValueError('Known implicit method required')
    stepper = classes[method](system, initial, initial_field,
        absolute_tolerance=spec.get('absolute_tolerance', 1e-12),
        relative_tolerance=spec.get('relative_tolerance', 1e-10))
    duration = (ref.period if spec['kind'] == 'contact' else 2*np.pi/abs(ref.omega))*spec.get('periods', 1.)
    count = spec.get('steps', 120); dt = duration/count
    state, field = initial.copy(), initial_field
    history = []; residuals = []; limits = []; iterations = []; evaluations = []; newton = []
    max_induction = 0.; max_refinement = 0
    for step in range(count):
        previous, previous_field = state.copy(), field
        try:
            state, field, diagnostic = stepper.step(state, field, dt)
        except Exception:
            if progress is not None:
                progress(dict(step=step, time=step*dt, failed_next_step=True,
                    density_min=float(state[..., 0].min()), max_residual=max(residuals, default=0.),
                    max_iterations=max(iterations, default=0)),
                    dict(state=state.copy(), magnetic=field.packed_increment(system.geometry)))
            raise
        residuals.append(diagnostic['residual']); iterations.append(diagnostic['iterations'])
        limits.append(diagnostic['residual_limit']); newton.append(diagnostic['newton_corrections'])
        evaluations.append(diagnostic['rhs_evaluations'])
        max_induction = max(max_induction, diagnostic['max_induction_residual'])
        max_refinement = max(max_refinement, diagnostic['max_induction_refinements'])
        if (step+1) % max(1, count//10) == 0 or step+1 == count:
            history.append(dict(step=step+1, time=(step+1)*dt, density_min=float(state[..., 0].min()),
                grid_change=grid_change(d, system, ref, initial, initial_field, state, field),
                max_residual=max(residuals), max_iterations=max(iterations), newton_corrections=sum(newton)))
            if progress is not None:
                progress(history[-1], dict(state=state.copy(), magnetic=field.packed_increment(system.geometry)))
    change = grid_change(d, system, ref, initial, initial_field, state, field)
    result = dict(spec=spec, status='complete', method='strong_wall_'+method, duration=duration, dt=dt,
        initial=field_norms(spec, ref, d, initial, system, initial_field, 0., 12),
        final=field_norms(spec, ref, d, state, system, field, duration, 12),
        grid_change=change, history=history, max_residual=max(residuals),
        max_iterations=max(iterations), rhs_evaluations=sum(evaluations),
        newton_corrections=sum(newton), max_residual_fraction=max(np.array(residuals)/np.array(limits)),
        max_induction_residual=max_induction, max_induction_refinements=max_refinement,
        coordinate_density=stepper.coordinates.density, coordinate_speed=stepper.coordinates.speed,
        elapsed_seconds=time.monotonic()-start)
    initial_energy = float(np.sum(d.volume*d.fields(initial)['energy']))
    model = RegularCurrent(system.mapping, initial[..., 0], ref.ion_mass)
    initial_energy += sum(initial_field.energy_components(model).values())
    initial_energy += float(np.sum(model.mesh.volumes*model.current_energy(initial_field.current(model))))
    result['grid_energy_relative_drift'] = change['total']/initial_energy
    if spec['kind'] == 'wave': result['energy_change_over_wave_energy'] = change['total']/ref.perturbation_energy()
    arrays = dict(initial_state=initial, final_state=state, penultimate_state=previous,
        initial_magnetic=initial_field.packed_increment(system.geometry),
        final_magnetic=field.packed_increment(system.geometry),
        penultimate_magnetic=previous_field.packed_increment(system.geometry),
        dt=np.array(dt), step_residuals=np.array(residuals), step_residual_limits=np.array(limits),
        step_iterations=np.array(iterations), step_newton_corrections=np.array(newton),
        step_evaluations=np.array(evaluations))
    return result, arrays
