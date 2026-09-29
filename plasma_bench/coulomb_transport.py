"""Conditional weakly coupled Maxwellian collision scales in SI units."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0, mu_0 as MU0, c

SOURCE="https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf"


def coulomb_coefficients(electron_density,electron_temperature_ev,ion_temperature_ev,ion_mass,charge_state=1):
    """NRL thermal Coulomb logarithms/rates; no clipping of invalid logarithms."""
    ne,te,ti,mi=map(float,(electron_density,electron_temperature_ev,ion_temperature_ev,ion_mass))
    if not all(np.isfinite(x) and x>0 for x in (ne,te,ti,mi)):
        raise ValueError("Finite positive charged-plasma density, temperatures and ion mass required")
    if isinstance(charge_state,bool) or not isinstance(charge_state,int) or charge_state<1:
        raise ValueError("Positive integer ion charge state required")
    z=charge_state;ni=ne/z
    if ti*ME/mi>=min(te,10*z*z):
        raise ValueError("Electron-ion Coulomb-log branch requires slow thermal ions")
    log_ei=(23-.5*np.log(ne/1e6)-np.log(z)+1.5*np.log(te) if te<10*z*z
            else 24-.5*np.log(ne/1e6)+np.log(te))
    log_ii=23-.5*np.log(2*ni/1e6)-3*np.log(z)+1.5*np.log(ti)
    if min(log_ei,log_ii)<=0:
        raise ValueError("Nonpositive Coulomb logarithm: classical closure unavailable, no floor applied")
    coupling=E**4/(4*np.pi*EPS0)**2
    nue=4*np.sqrt(2*np.pi)*ne*z*coupling*log_ei/(3*np.sqrt(ME)*(E*te)**1.5)
    nui=4*np.sqrt(np.pi)*ni*z**4*coupling*log_ii/(3*np.sqrt(mi)*(E*ti)**1.5)
    eta=ME*nue/(ne*E**2)
    return {"coulomb_log_ei":float(log_ei),"coulomb_log_ii":float(log_ii),
            "electron_ion_momentum_rate_s":float(nue),"ion_ion_collision_rate_s":float(nui),
            "momentum_resistivity_ohm_m":float(eta),
            "parallel_resistivity_z1_ohm_m":float(eta/1.96) if z==1 else None,
            "electron_thermal_mean_free_path_m":float(np.sqrt(E*te/ME)/nue),
            "ion_thermal_mean_free_path_m":float(np.sqrt(E*ti/mi)/nui)}


def regime(electron_density,electron_temperature_ev,ion_temperature_ev,ion_mass,*,
           charge_state=1,ionization=1.,field=0.,current_density=0.,length=1.,timescale=1.):
    """Local/reference diagnostics, never a sufficient declaration of validity.

    length is the smallest physical structure/gradient scale being assessed,
    not automatically the grid spacing. timescale is a physical variation time.
    """
    ne,te,ti,mi,b,j,l,t=map(float,(electron_density,electron_temperature_ev,ion_temperature_ev,
                                 ion_mass,abs(field),abs(current_density),length,timescale))
    if not all(np.isfinite(x) for x in (ne,te,ti,mi,b,j,l,t,ionization)) or min(ne,te,ti,mi,l,t)<=0:
        raise ValueError("Finite positive plasma state, length and timescale required")
    if isinstance(charge_state,bool) or not isinstance(charge_state,int) or charge_state<1 or not 0<ionization<=1:
        raise ValueError("Positive integer charge and ionization in (0,1] required")
    z=charge_state;ni=ne/z;ve=np.sqrt(E*te/ME);vi=np.sqrt(E*ti/mi)
    de=np.sqrt(ME/(MU0*ne*E**2));di=np.sqrt(mi/(MU0*ni*(z*E)**2))
    ld=np.sqrt(EPS0*E*te/(ne*E**2));drift=j/(ne*E);omega=E*b/ME
    warnings=[]
    out={"source":SOURCE,"scope":"Conditional local Maxwellian Coulomb estimates; no total-neutral or kinetic closure",
         "electron_density_m3":ne,"electron_temperature_ev":te,"ion_temperature_ev":ti,
         "ion_mass_kg":mi,"charge_state":z,"ionization_fraction":float(ionization),
         "field_t":b,"current_density_a_m2":j,"structure_length_m":l,"variation_time_s":t,
         "electron_inertial_length_m":float(de),"ion_inertial_length_m":float(di),
         "electron_debye_length_m":float(ld),"electron_inertia_over_length":float(de/l),
         "hall_length_ratio":float(di/l),"debye_over_length":float(ld/l),
         "electron_plasma_period_scale_s":float(de/c),"light_transit_over_time":float(l/(c*t)),
         "current_drift_m_s":float(drift),"drift_over_electron_thermal_speed":float(drift/ve),
         "drift_over_ion_thermal_speed":float(drift/vi),"drift_over_light_speed":float(drift/c),
         "electron_gyro_radius_m":float(ve/omega) if omega else None,
         "classical_coefficients":None}
    for ratio,label in ((di/l,"Hall scale"),(de/l,"electron inertia scale"),(ld/l,"Debye scale"),
                        (drift/ve,"electron drift/thermal ratio"),(drift/vi,"ion-relative drift/thermal ratio")):
        if ratio>=.1:warnings.append(f"{label} is not small (>=0.1); review the reduced model locally.")
    if drift/c>=.1:warnings.append("Current implies relativistic carrier drift; this nonrelativistic model is inadequate.")
    if l/(c*t)>=.01:warnings.append("Electromagnetic transit is not negligible; quasistatic vacuum needs review.")
    if ionization<1:warnings.append("Neutral collisions, ionization and recombination are absent; Coulomb rates are only one contribution.")
    try:
        coefficients=coulomb_coefficients(ne,te,ti,mi,z);out["classical_coefficients"]=coefficients
        nue=coefficients["electron_ion_momentum_rate_s"];nui=coefficients["ion_ion_collision_rate_s"]
        out.update(electron_collision_time_over_variation=float(1/(nue*t)),
                   ion_collision_time_over_variation=float(1/(nui*t)),electron_hall_parameter=float(omega/nue),
                   electron_mean_free_path_over_length=coefficients["electron_thermal_mean_free_path_m"]/l,
                   ion_mean_free_path_over_length=coefficients["ion_thermal_mean_free_path_m"]/l)
        if min(coefficients["coulomb_log_ei"],coefficients["coulomb_log_ii"])<5:
            warnings.append("Small Coulomb logarithm (<5); weak-coupling transport is poorly controlled.")
        if max(out["electron_collision_time_over_variation"],out["ion_collision_time_over_variation"])>=.1:
            warnings.append("Thermal collisions are slow relative to variation; local collisional closure needs review.")
        if max(out["electron_mean_free_path_over_length"],out["ion_mean_free_path_over_length"])>=.1:
            warnings.append("Thermal mean free path is not small; anisotropic/nonlocal or kinetic transport may be needed.")
        if omega/nue>=1:warnings.append("Electrons are magnetized; scalar momentum resistivity is not a full conductivity tensor.")
        if omega and ve/omega<=10*ld:warnings.append("Electron gyro radius is not well above Debye length; classical collision assumptions need review.")
    except ValueError as error:
        warnings.append(str(error))
    out["warnings"]=warnings
    out["validity_established"]=False
    return out
