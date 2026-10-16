"""One collision force/heat pair for the moving finite-mass plasma."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME
from .coulomb_exchange_wall_plasma import equilibration_rate
from .maxwellian_pair_collisions import pair_moments


def collision_config(model, fluid, resistivity, coulomb_exchange):
    """Validate the optional paired closure without adding parallel collisions."""
    if model is None:
        return None
    if (not isinstance(model, dict) or set(model) != {'kind', 'coulomb_log'}
            or model['kind'] != 'maxwellian_pair'):
        raise ValueError('Collision model requires kind=maxwellian_pair and explicit coulomb_log')
    log = model['coulomb_log']
    if not np.isscalar(log) or not np.isfinite(log) or log <= 0:
        raise ValueError('A positive finite declared Coulomb logarithm is required')
    if resistivity != 0 or coulomb_exchange:
        raise ValueError('Paired collisions replace prescribed resistivity and separate thermal exchange')
    if not np.all(fluid.gammas == 5/3):
        raise ValueError('Isotropic Maxwellian collisions require monatomic gamma=5/3 species')
    return dict(kind='maxwellian_pair', coulomb_log=float(log))


def collision_sources(fields, drift, ion_mass, *, model, resistivity, coulomb_exchange):
    """Return electric forcing R_e/(e*n) and the same pair's species heat.

    The moving solver uses w=Ui-Ue and J=e*n*w. Temperatures and density
    belong to positive plasma quadrature; this function never fills vacuum.
    """
    n = fields['n'].ravel(); te, ti = fields['temperature'].reshape(-1, 2).T
    if model is None:
        electric = resistivity*E*n[:, None]*drift
        drag_power = resistivity*(E*n)**2*np.sum(drift*drift, axis=-1)
        frequency = np.zeros_like(n); exchange = np.zeros_like(n); log = None
        if coulomb_exchange:
            frequency, log = equilibration_rate(n, te, ti, ion_mass)
            exchange = 1.5*E*n*frequency*(ti-te)
        share = ion_mass/(ion_mass+ME)
        heating = np.column_stack((share*drag_power+exchange, (1-share)*drag_power-exchange))
        relaxation = resistivity*E**2*n*(1/ME+1/ion_mass)
    else:
        pair = pair_moments(n, n, te, ti, -drift, ME, ion_mass,
                            coulomb_log=model['coulomb_log'])
        electric = pair['force_on_a']/(E*n[:, None])
        heating = np.column_stack((pair['heat_to_a'], pair['heat_to_b']))
        drag_power = pair['drag_power']; exchange = pair['thermal_exchange_to_a']
        relaxation = pair['drag_coefficient']/n*(1/ME+1/ion_mass)
        frequency = None; log = np.full_like(n, model['coulomb_log'])
    return dict(electric=electric, heating=heating, drag_power=drag_power,
                thermal_exchange_to_electrons=exchange, relaxation_rate=relaxation,
                exchange_frequency=frequency, coulomb_log=log)
