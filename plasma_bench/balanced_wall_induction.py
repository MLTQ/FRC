"""Strict kernel-separated induction using the declared wall Galerkin curl."""
from copy import copy
import numpy as np
from .inertial_induction import InductionRate
from .balanced_wall_core import BalancedWallCore
from .uniform_wall_plasma import UniformWallInduction
from .kernel_induction import separate, difference_curl


class BalancedWallInduction(UniformWallInduction):
    def solve(self, forcing):
        model = self.model; geo = model.geometry; f = np.asarray(forcing, float)
        if f.shape != (model.mass.shape[0],) or not np.isfinite(f).all():
            raise ValueError('Finite complete weak forcing required')
        electric0 = self.mass_inverse(f); small, kernel = separate(geo, electric0)
        weak = model.mass@small; norm = np.linalg.norm; attempts = []; self.attempt_trace = []
        relative = lambda a, b: float(norm(a)/max(norm(b), np.finfo(float).tiny))
        for refinement in range(5):
            inner = copy(self); inner.tolerance = self.tolerance*10.**(-refinement)
            try: result = BalancedWallCore.solve(inner, weak)
            except RuntimeError as exc:
                self.attempt_trace.append(dict(tolerance=inner.tolerance, trace=getattr(inner, 'refinement_trace', []), error=str(exc)))
                attempts.append(dict(inner_tolerance=inner.tolerance, error=str(exc))); continue
            self.attempt_trace.append(dict(tolerance=inner.tolerance, trace=inner.refinement_trace))
            electric = result.electric+kernel; inertia = model.inertia@result.current
            original_curl = difference_curl(geo, electric); curl = geo.project_curl(original_curl)
            ohm = relative(model.mass@electric-inertia-f, abs(model.mass@electric)+abs(inertia)+abs(f))
            faraday = relative(curl-result.magnetic, abs(geo.faraday)@abs(electric)+abs(result.magnetic))
            reconstruction = relative(electric0-small-kernel, abs(electric0)+abs(small)+abs(kernel))
            components = dict(result.residual_components, full_ohm_backward=ohm,
                full_faraday_backward=faraday, kernel_reconstruction_backward=reconstruction)
            residual = max(components.values())
            attempts.append(dict(inner_tolerance=inner.tolerance, residual=residual))
            if np.isfinite(residual) and residual <= self.tolerance:
                components.update(inner_tolerance=inner.tolerance, outer_refinements=refinement,
                    full_faraday_rate_relative=relative(curl-result.magnetic, result.magnetic),
                    unprojected_curl_relative=relative(original_curl-result.magnetic, result.magnetic),
                    kernel_to_remainder_ratio=float(norm(kernel)/norm(small)) if norm(small)>0 else None)
                return InductionRate(result.magnetic, result.current, electric, result.iterations,
                    residual, result.refinements, components)
        raise RuntimeError(f'Strong-wall induction residual failure: {attempts}; no input modified')
