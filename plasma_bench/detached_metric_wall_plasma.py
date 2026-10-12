"""Geometry cache templates that do not retain a completed density model."""
import gc
from .direct_metric_wall_plasma import DirectMetricWallPlasma


class DetachedMetricWallPlasma(DirectMetricWallPlasma):
    def operators(self, model):
        result = super().operators(model)
        if self._operators is not None:
            # Returned shallow copies already reference the actual current model.
            # Templates retain maps/factors, not a completed state's density.
            for template in self._operators:
                template.model = None
        return result

    def derivative(self, state, field):
        result = super().derivative(state, field)
        # LinearOperator bound methods form a model/inertia reference cycle.
        # Collect only after the complete parent derivative has returned.
        gc.collect()
        return result


def detached_metric(system):
    result = DetachedMetricWallPlasma(system.mapping, system.ion_mass,
        cache_geometry=system.cache_geometry, **system.solver_options)
    result.original_system = system.original_system
    return result
