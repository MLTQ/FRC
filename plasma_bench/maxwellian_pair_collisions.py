"""Paired Landau drag and heat for two displaced isotropic Maxwellians."""
import numpy as np
from scipy.constants import elementary_charge as E, epsilon_0 as EPS0
from scipy.special import erf, gammainc


def pair_moments(density_a, density_b, temperature_a_ev, temperature_b_ev,
                 drift_a_minus_b, mass_a, mass_b, *, coulomb_log,
                 charge_a=-1, charge_b=1):
    """Return force on species a [N/m³] and species heat [W/m³].

    Temperature is in eV, drift in m/s, density in m^-3, mass in kg.
    The positive supplied logarithm is a declared weak-scattering approximation;
    this function does not infer a finite-drift screening/cutoff prescription.
    Exact instantaneous Maxwellian moments do not preserve Maxwellian shapes.
    """
    u = np.asarray(drift_a_minus_b, float)
    if u.ndim < 1 or u.shape[-1] != 3 or not np.isfinite(u).all():
        raise ValueError('Finite three-component relative species velocity required')
    values = np.broadcast_arrays(density_a, density_b, temperature_a_ev,
                                temperature_b_ev, coulomb_log, u[..., 0])
    na, nb, tea, teb, log, _ = [np.asarray(v, float) for v in values]
    if (not all(np.isfinite(v).all() and np.all(v > 0) for v in (na, nb, tea, teb, log))
            or not np.isfinite([mass_a, mass_b]).all() or min(mass_a, mass_b) <= 0):
        raise ValueError('Positive finite charged-plasma states, masses and Coulomb logarithm required')
    if any(isinstance(z, bool) or not isinstance(z, int) or z == 0 for z in (charge_a, charge_b)):
        raise ValueError('Nonzero integer charge numbers required')
    u = np.broadcast_to(u, na.shape+(3,))
    ta, tb = E*tea, E*teb
    total_mass = mass_a+mass_b; reduced_mass = mass_a*mass_b/total_mass
    variance = ta/mass_a+tb/mass_b; speed2 = np.sum(u*u, axis=-1)
    x2 = speed2/(2*variance); x = np.sqrt(x2)
    # Analytic zero-drift limits, with Taylor evaluation only to avoid roundoff
    # cancellation. These are not physical floors or a softened collision law.
    small = x < .05
    a_small = np.sqrt(2/np.pi)/(3*variance**1.5)*(1+x2*(-3/5+x2*(3/14+x2*(-1/18+x2/88))))
    b_small = np.sqrt(2/np.pi)/np.sqrt(variance)*(1+x2*(-1/3+x2*(1/10+x2*(-1/42+x2/216))))
    a_large = np.divide(gammainc(1.5, x2), speed2**1.5,
                        out=np.zeros_like(speed2), where=~small)
    b_large = np.divide(erf(x), np.sqrt(speed2), out=np.zeros_like(speed2), where=~small)
    coefficient_a = np.where(small, a_small, a_large)
    potential_b = np.where(small, b_small, b_large)
    coupling = (charge_a*charge_b*E**2/(4*np.pi*EPS0))**2
    prefactor = 4*np.pi*na*nb*coupling*log/reduced_mass
    drag_coefficient = prefactor*coefficient_a
    force = -drag_coefficient[..., None]*u
    drag_power = drag_coefficient*speed2
    exchange = prefactor*(tb-ta)/(total_mass*variance)*potential_b
    qa = (ta/mass_a)/variance*drag_power+exchange
    qb = (tb/mass_b)/variance*drag_power-exchange
    if not all(np.isfinite(v).all() for v in (force, qa, qb, drag_coefficient)):
        raise ValueError('Nonfinite collision moment; no state repair')
    return dict(force_on_a=force, force_on_b=-force, heat_to_a=qa, heat_to_b=qb,
                drag_power=drag_power, drag_coefficient=drag_coefficient,
                thermal_exchange_to_a=exchange,
                entropy_production_over_kB=qa/ta+qb/tb,
                relative_drift_over_combined_thermal_speed=np.sqrt(speed2/variance))
