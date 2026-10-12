"""Record, without repair, the arithmetic paths of one strict iterate."""
import numpy as np


def observe(solver, solution, e0, rate, current, electric, rhs, operator_value, difference):
    geo = solver.model.geometry; inertia = solver.model.inertia
    D = solver.dual_lift; T = solver.lift; F = geo.faraday
    norm = np.linalg.norm
    scale = max(norm(rhs), np.finfo(float).tiny)
    pre = D@solution; post = geo.hessian@rate
    jpre = solver.mass_inverse(F.T@pre)
    epre = solver.mass_inverse(inertia@jpre)
    ei = solver.mass_inverse(inertia@current)
    f0 = F@e0; fi = F@ei; fp = F@epre; full = F@electric
    metric = T.T@pre; projected_rate = D.T@rate
    op = operator_value-rhs; dual = D.T@difference; gap = op+dual
    electric_gap = epre+ei
    z = electric-e0
    addition_low = (e0-(electric-z))+(ei-z)
    terms = dict(
        metric_association=metric-projected_rate,
        inertia_paths=D.T@(F@electric_gap),
        inertia_curl_distribution=D.T@((fp+fi)-F@electric_gap),
        inertia_dual_distribution=(D.T@fp+D.T@fi)-D.T@(fp+fi),
        electric_curl_distribution=D.T@(full-f0-fi),
        electric_dual_distribution=(D.T@full-rhs-D.T@fi)-D.T@(full-f0-fi),
        rate_dual_distribution=dual-(D.T@full-projected_rate))
    predicted = sum(terms.values()); remainder = gap-predicted
    addition_prediction = -D.T@(F@addition_low)
    arrays = dict(solution=solution, rhs=rhs, magnetic_rate=rate, current_rate=current,
        forcing_electric=e0, inertia_electric=ei, electric=electric,
        operator_residual=op, physical_dual_residual=dual, opposing_gap=gap,
        predicted_gap=predicted, decomposition_remainder=remainder,
        dual_action_gap=pre-post, current_path_gap=jpre+current,
        inertia_electric_gap=electric_gap, electric_addition_low=addition_low,
        electric_addition_prediction=addition_prediction, **terms)
    summary = dict(opposing_gap_over_rhs=float(norm(gap)/scale),
        terms_over_rhs={k: float(norm(v)/scale) for k, v in terms.items()},
        decomposition_remainder_over_rhs=float(norm(remainder)/scale),
        decomposition_remainder_fraction=float(norm(remainder)/max(norm(gap), np.finfo(float).tiny)),
        electric_addition_prediction_over_rhs=float(norm(addition_prediction)/scale),
        curl_distribution_minus_addition_fraction=float(norm(terms['electric_curl_distribution']-addition_prediction)/
            max(norm(terms['electric_curl_distribution']), np.finfo(float).tiny)),
        dual_action_relative=float(norm(pre-post)/max(norm(post), np.finfo(float).tiny)),
        current_action_relative=float(norm(jpre+current)/max(norm(current), np.finfo(float).tiny)))
    return summary, {k: v.copy() for k, v in arrays.items()}
