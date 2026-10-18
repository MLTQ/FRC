"""First-cut field-reversed theta-pinch design relations: explicit scaling, not a simulation."""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.constants import mu_0, k as KB, e as QE, epsilon_0, atomic_mass as AMU, c as C0, torr as TORR

# Atomic mass [u] of the dissociated atom and the ion charge assumed after formation.
GASES = {"H2": (1.00784, 1), "D2": (2.01410, 1), "N2": (14.0067, 1)}


def rigid_rotor_flux(k):
    """phi_p/(pi r_s^2 B_e) for Bz = B_e tanh(K(2r^2/r_s^2-1)); null at r_s/sqrt2."""
    return np.log(np.cosh(k))/(2*k)


def rigid_rotor_decay(k):
    """tau_phi*eta/(mu_0 r_s^2): flux over Ohmic loss rate at the null, fixed profile."""
    return np.log(np.cosh(k))/(8*k*k)


def coulomb_log(ne_m3, te_ev, z=1):
    """NRL electron-ion Coulomb logarithm."""
    n_cc = ne_m3*1e-6
    return np.where(te_ev < 10*z*z, 23-np.log(np.sqrt(n_cc)*z*te_ev**-1.5),
                    24-np.log(np.sqrt(n_cc)/te_ev))


def spitzer_perp(te_ev, ln_lambda, z=1):
    """NRL perpendicular Spitzer resistivity [ohm m]; about twice the parallel value."""
    return 1.03e-4*z*ln_lambda*te_ev**-1.5


def sheet_coil_center_field(current, radius, length):
    """Axial field at the centre of a thin single-layer solenoid carrying `current` ampere-turns."""
    return mu_0*current/np.hypot(length, 2*radius)


def wheeler_inductance(radius, length, turns=1):
    """Wheeler's short-solenoid inductance, mu_0 pi a^2 N^2/(l+0.9a); about 1% for l > 0.8a."""
    return mu_0*np.pi*radius**2*turns**2/(length+.9*radius)


@dataclass(frozen=True)
class Assumptions:
    """Empirical formation efficiencies: round mid-range values, not fitted to any device."""
    particle_retention: float = .4   # trapped ions / fill atoms in the coil bore
    flux_retention: float = .4       # trapped poloidal flux / initial bias flux in the bore
    separatrix_to_wall: float = .5   # r_s / r_w
    length_to_coil: float = .6       # l_s / l_c
    shape: float = .75               # V_s / (pi r_s^2 l_s); ellipsoid 2/3, racetrack ~1
    profile_k: float = 2.            # rigid-rotor steepness; edge/peak density sech^2 K
    ti_over_te: float = 1.


@dataclass(frozen=True)
class Device:
    wall_radius: float = .071        # Pyrex bore [m]
    coil_radius: float = .080        # single-turn sheet coil [m]
    coil_length: float = .5          # [m]
    turns: int = 1                   # series turns per segment, wound uniformly
    segments: int = 1                # parallel coil sections, each fed by its own module
    stray_inductance: float = 30e-9  # bank, switch and feed, per module [H]
    gas_temperature_k: float = 300.
    rise_per_implosion: float = 3.   # main-bank quarter period / snowplow time
    pi_skin_margin: float = 10.      # required PI column skin time / quarter period


def pi_temperature(n_e, skin_time, wall_radius, z=1):
    """Electron temperature [eV] whose column skin time mu_0 r_w^2/eta equals `skin_time`."""
    te = 2.
    for _ in range(50):
        te = (1.03e-4*z*coulomb_log(n_e, te, z)*skin_time/(mu_0*wall_radius**2))**(2/3)
    return float(te)


def design_point(gas, fill_mtorr, b_ext, device=Device(), a=Assumptions()):
    """Equilibrium, formation and drive estimates for one gas/fill/external-field choice."""
    mass_u, z = GASES[gas]
    m_i = mass_u*AMU
    r_w, r_c, l_c = device.wall_radius, device.coil_radius, device.coil_length
    n_fill = 2*fill_mtorr*1e-3*TORR/(KB*device.gas_temperature_k)   # atoms, diatomic fill
    r_s = a.separatrix_to_wall*r_w
    x_s = r_s/r_c
    beta = 1-x_s**2/2                                              # elongated, flux-conserving
    volume = a.shape*np.pi*r_s**2*a.length_to_coil*l_c
    n_avg = a.particle_retention*n_fill*np.pi*r_w**2*l_c/volume
    n_peak = n_avg/beta                                            # isothermal profile
    p_mag = b_ext**2/(2*mu_0)
    te = beta*p_mag/(n_avg*QE*(a.ti_over_te+z))                     # n_i(Ti+Z Te) = <beta> p_mag
    ti = a.ti_over_te*te
    ln_l = coulomb_log(z*n_peak, te, z)
    eta = spitzer_perp(te, ln_l, z)
    rho_i = np.sqrt(2*m_i*QE*ti)/(z*QE*b_ext)
    d_i = C0/np.sqrt(n_peak*(z*QE)**2/(epsilon_0*m_i))
    phi_p = rigid_rotor_flux(a.profile_k)*np.pi*r_s**2*b_ext
    v_a = b_ext/np.sqrt(mu_0*n_fill*m_i)
    t_sp = np.sqrt(2)*r_w/v_a                                      # planar snowplow transit
    current = b_ext*np.hypot(l_c, 2*r_c)/mu_0                      # total ampere-turns
    l_coil = wheeler_inductance(r_c, l_c)                          # single-turn equivalent
    k, n = device.segments, device.turns
    # Synchronous segments share the coil flux: each module sees k N^2 L_coil plus its own stray.
    l_mod = k*n*n*l_coil+device.stray_inductance
    i_mod = current/(k*n)
    t_q = device.rise_per_implosion*t_sp
    omega = np.pi/(2*t_q)
    return dict(
        gas=gas, fill_mtorr=fill_mtorr, b_ext_t=b_ext, fill_atoms_m3=n_fill,
        separatrix_radius_m=r_s, x_s=x_s, avg_beta=beta, n_avg_m3=n_avg, n_peak_m3=n_peak,
        compression=n_avg/n_fill, te_ev=te, ti_ev=ti, total_t_ev=te+ti, coulomb_log=ln_l,
        eta_ohm_m=eta, tau_phi_classical_s=mu_0*r_s**2*rigid_rotor_decay(a.profile_k)/eta,
        rho_i_m=rho_i, s_bar=r_s/(4*rho_i), rs_over_di=r_s/d_i, trapped_flux_wb=phi_p,
        bias_t=-phi_p/(a.flux_retention*np.pi*r_w**2), snowplow_s=t_sp,
        quarter_period_s=t_q, coil_current_a=current, coil_inductance_h=l_coil,
        turns=n, segments=k, module_current_a=i_mod, module_inductance_h=l_mod,
        loop_voltage_v=l_coil*current*omega, bank_energy_j=k*l_mod*i_mod**2/2,
        bank_voltage_v=i_mod*l_mod*omega, module_capacitance_f=1/(omega**2*l_mod),
        bank_capacitance_f=k/(omega**2*l_mod),
        pi_te_min_ev=pi_temperature(z*n_fill, device.pi_skin_margin*t_q, r_w, z))


def record(device, assumptions):
    return dict(device=asdict(device), assumptions=asdict(assumptions))
