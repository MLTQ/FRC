"""Strict inner refinement controlled by every original full induction residual."""
from copy import copy
import numpy as np
from .inertial_induction import InertialInduction, InductionRate
from .kernel_induction import KernelInduction, separate, difference_curl


class RefinedKernelInduction(KernelInduction):
    """Same equations and outer tolerance; tighten only the inner solve as needed."""
    def solve(self, forcing):
        model = self.model; f = np.asarray(forcing, float)
        if f.shape != (model.mass.shape[0],) or not np.isfinite(f).all():
            raise ValueError('Finite full weak electric forcing required')
        electric0 = self.mass_inverse(f)
        small, kernel = separate(model.geometry, electric0)
        weak = model.mass@small
        norm = np.linalg.norm
        relative = lambda a, b: float(norm(a)/max(norm(b), np.finfo(float).tiny))
        attempts = []
        for refinement in range(5):
            inner = copy(self); inner.tolerance = self.tolerance*10.**(-refinement)
            try:
                result = InertialInduction.solve(inner, weak)
            except RuntimeError as exc:
                attempts.append(dict(inner_tolerance=inner.tolerance, error=str(exc)))
                continue
            electric = result.electric+kernel
            inertia = model.inertia@result.current
            ohm = relative(model.mass@electric-inertia-f, abs(model.mass@electric)+abs(inertia)+abs(f))
            curl = difference_curl(model.geometry, electric)
            bound = abs(model.geometry.faraday)@abs(electric)+abs(result.magnetic)
            faraday = relative(curl-result.magnetic, bound)
            reconstruction = relative(electric0-small-kernel, abs(electric0)+abs(small)+abs(kernel))
            components = dict(result.residual_components, full_ohm_backward=ohm,
                full_faraday_backward=faraday, kernel_reconstruction_backward=reconstruction)
            residual = max(components.values())
            attempts.append(dict(inner_tolerance=inner.tolerance, residual=residual,
                                 full_faraday_backward=faraday))
            if np.isfinite(residual) and residual <= self.tolerance:
                components['full_faraday_rate_relative'] = relative(curl-result.magnetic, result.magnetic)
                components['kernel_to_remainder_ratio'] = float(norm(kernel)/norm(small)) if norm(small)>0 else None
                components['inner_tolerance'] = inner.tolerance
                components['outer_refinements'] = refinement
                return InductionRate(result.magnetic, result.current, electric, result.iterations,
                                     residual, result.refinements, components)
        raise RuntimeError(f'All original induction residuals failed after strict inner refinement: {attempts}; no input modified')
