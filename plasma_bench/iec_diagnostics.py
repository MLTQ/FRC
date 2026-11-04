"""Full-population observables; no symmetry restriction on the evolved plasma."""
import numpy as np
from scipy.constants import elementary_charge as qe


def ion_motion(p):
    select = p.species == 1
    x, v = p.x[select], p.v[select]
    count = len(x)
    if not count:
        return dict(count=0, radius_rms_m=None, radial_speed_m_s=None)
    centre, velocity = x.mean(axis=0), v.mean(axis=0)
    dx, dv = x-centre, v-velocity
    r2 = float(np.mean(np.sum(dx*dx, axis=1)))
    radius = np.sqrt(r2)
    correlation = float(np.mean(np.sum(dx*dv, axis=1)))
    speed = correlation/radius if radius > 0 else 0.
    alpha = correlation/r2 if r2 > 0 else 0.
    factor = .5*p.ion_mass*p.weight
    translation = factor*count*float(velocity@velocity)
    breathing = factor*alpha**2*float(np.sum(dx*dx))
    residual = factor*float(np.sum((dv-alpha*dx)**2))
    total = factor*float(np.sum(v*v))
    principal = np.sqrt(np.maximum(np.linalg.eigvalsh(dx.T@dx/count), 0.))
    return dict(count=count, radius_rms_m=float(radius), radial_speed_m_s=float(speed),
                centroid_m=centre.tolist(), principal_rms_m=principal.tolist(),
                radius_origin_rms_m=float(np.sqrt(np.mean(np.sum(x*x, axis=1)))),
                radius_quantiles_m=np.quantile(np.linalg.norm(dx, axis=1), [.1,.5,.9]).tolist(),
                translation_j=translation, breathing_j=breathing, residual_motion_j=residual,
                kinetic_decomposition_error_j=translation+breathing+residual-total,
                mean_energy_ev=total/(count*p.weight*qe))


def ion_sensor(p):
    """Cheap per-step radius and dilation speed; diagnostics use the same definition."""
    x, v = p.x[p.species == 1], p.v[p.species == 1]
    if not len(x):
        return None, 0., 0
    dx, dv = x-x.mean(axis=0), v-v.mean(axis=0)
    radius = np.sqrt(np.mean(np.sum(dx*dx, axis=1)))
    return float(radius), float(np.mean(np.sum(dx*dv, axis=1))/radius) if radius > 0 else 0., len(x)
