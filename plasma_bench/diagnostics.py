"""Reference scales and candid model-validity warnings."""
import numpy as np
from scipy.constants import mu_0, elementary_charge as E, electron_mass as ME, c
from .coulomb_transport import regime


def characteristic_scales(config):
    m, g = config.material, config.grid
    b, length, sigma = abs(config.bias), m.region.radius, m.conductivity
    ni, ne = m.ion_density, m.ion_density*m.charge_state
    pressure = E*(ni*m.ion_temperature_ev+ne*m.electron_temperature_ev)
    va = b/np.sqrt(mu_0*ni*m.ion_mass)
    di = np.sqrt(m.ion_mass/(mu_0*ni*(m.charge_state*E)**2))
    wi, we = m.charge_state*E*b/m.ion_mass, E*b/ME
    vi, ve = np.sqrt(2*E*m.ion_temperature_ev/m.ion_mass), np.sqrt(2*E*m.electron_temperature_ev/ME)
    nu_eff = ne*E**2/(ME*sigma)
    out = dict(reference_field_t=b, structure_length_m=length, pressure_pa=pressure,
               beta=None if b == 0 else 2*mu_0*pressure/b**2,
               alfven_speed_m_s=va, alfven_time_s=None if va == 0 else length/va,
               magnetic_reynolds_number=0.0, lundquist_number=mu_0*sigma*length*va,
               diffusion_time_s=mu_0*sigma*length**2,
               ion_skin_depth_m=di, ion_skin_depth_over_size=di/length,
               ion_gyrofrequency_rad_s=wi, electron_gyrofrequency_rad_s=we,
               ion_gyroradius_m=None if wi == 0 else vi/wi,
               electron_gyroradius_m=None if we == 0 else ve/we,
               inferred_drude_collision_rate_s=nu_eff,
               inferred_electron_hall_parameter=we/nu_eff,
               physical_collision_frequency_s=None, physical_mean_free_path_m=None,
               dr_m=g.dr, dz_m=g.dz, dr_over_ion_skin_depth=g.dr/di,
               plasma_radial_cells=length/g.dr, light_crossing_time_s=2*g.half_length/c)
    warnings = []
    if config.stage == 0:
        warnings.append("Stage 0 is vacuum-only; plasma scales below are reference inputs, not evolved plasma.")
    else:
        warnings.append("Stage 1 holds density, temperature and plasma position fixed; J×B motion and thermal feedback are omitted.")
        if out["alfven_time_s"] is not None and config.duration > out["alfven_time_s"]:
            warnings.append("Run exceeds the reference Alfven time; stationary plasma may be a poor approximation.")
    if di/length > .1:
        warnings.append("Ion skin depth exceeds 0.1 of plasma radius; single-fluid/Hall-free modeling is questionable.")
    if wi and vi/wi > .1*length:
        warnings.append("Ion gyroradius is not small compared with plasma radius; magnetized-fluid assumptions need review.")
    if we/nu_eff > 1:
        warnings.append("Inferred electron Hall parameter exceeds one; scalar conductivity needs a transport justification.")
    if m.initial_ionization < .9:
        warnings.append("Partially ionized input: neutral collisions and ionization dynamics are omitted.")
    warnings.extend(["Conductivity remains prescribed; conditional Coulomb estimates below are not a complete transport closure.",
                     "Scales use |bias field|; local validity can be worse near magnetic nulls.",
                     "No equilibrium, stability, FRC success or particle confinement is established by these stages."])
    if config.stage == 1:
        warnings.append("Induced flux is fixed to zero at the outer boundary; enlarge the domain before trusting boundary-sensitive results.")
    if min(length/g.dr, m.region.half_length/g.dz) < 10:
        warnings.append("Plasma has fewer than ten cells along a dimension; spatial convergence is essential.")
    if not b:
        warnings.append("Zero bias: Bz/Bbias and reversal relative to bias are undefined.")
    features, frequencies = [], []
    for coil in config.coils:
        wave = coil.waveform
        for attr in ("rise", "fall", "width", "decay"):
            if hasattr(wave, attr):
                features.append(getattr(wave, attr))
        if hasattr(wave, "frequency"):
            frequencies.append(wave.frequency)
            features.append(1/wave.frequency)
        if hasattr(wave, "times"):
            features.extend(np.diff(wave.times))
    if features:
        out["shortest_waveform_scale_s"] = float(min(features))
        out["dt_over_waveform_scale"] = float(np.diff(config.times).max()/min(features))
        out["light_crossing_over_waveform_scale"] = float(out["light_crossing_time_s"]/min(features))
        if out["dt_over_waveform_scale"] > .05*(1+1e-12):
            warnings.append("Fewer than 20 time steps resolve the shortest waveform scale; refine dt.")
        if out["light_crossing_over_waveform_scale"] > .01:
            warnings.append("Waveform is fast compared with light transit; magnetoquasistatics need review.")
    if frequencies:
        skin = np.sqrt(1/(np.pi*mu_0*sigma*max(frequencies)))
        out["oscillatory_skin_depth_m"] = float(skin)
        out["radial_cells_per_skin_depth"] = float(skin/g.dr)
        if skin/g.dr < 5:
            warnings.append("Oscillatory skin depth spans fewer than five radial cells; refine the mesh.")
    reference=regime(ne,m.electron_temperature_ev,m.ion_temperature_ev,m.ion_mass,
                     charge_state=m.charge_state,ionization=m.initial_ionization,field=b,
                     length=length,timescale=min(features) if features else config.duration)
    reference["current_assumption"]="Zero-current initial reference; evolving local current must be assessed separately"
    out["coulomb_reference"]=reference
    warnings.extend("Coulomb reference: "+w for w in reference["warnings"])
    return out, warnings


def magnetic_energy(grid, br, bz, region):
    r, z = grid.mesh
    density = (br*br+bz*bz)/(2*mu_0)*2*np.pi*r
    return float(np.trapezoid(np.trapezoid(np.where(region.contains(r, z), density, 0), grid.z, axis=1), grid.r))


def longest_sampled_interval(times, flags):
    start, longest = None, 0.0
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        if start is not None and (not flag or i == len(flags)-1):
            end = i if flag else i-1
            longest = max(longest, float(times[end]-times[start]))
            start = None
    return longest
