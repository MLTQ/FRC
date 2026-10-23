"""Static neutral reservoir ionized by plasma electrons, for partially pre-ionized fills (frc-ibl.4; see .md).

A Lie-split source after each ConservativeVacuumMHD step, built like density_floor.apply_floor:
    dN = n_n (1 - exp(-n_e K_iz(T) dt))           per cell (exact for n_e and T frozen over the step)
New ions are born at rest: mass rises by m_i dN, momentum is unchanged, so the mass-loading kinetic-energy drop
stays in the cell as heat (total energy unchanged by it). The ionization cost E_c per new ion (declared) is removed
from the cell's thermal energy, capped so at most half of it is spent in one step. T is the single-fluid
temperature (p = 2 n T). The neutrals do not move, drag, or charge-exchange.
"""
import numpy as np
from scipy.constants import e as QE
from .conservative_coupled import ConservativeVacuumMHD, ConservativeState
from .molecular_beb import MolecularBEB


def beb_rate(directory='data/molecular/nist_beb', gas='H2', t_min=.3, t_max=2000., points=60):
    """Maxwellian H2 total-ionization rate K(T) [m^3/s] from NIST BEB, interpolated in log-log; 0 below t_min."""
    beb = MolecularBEB.from_nist(directory, gas)
    t = np.geomspace(t_min, t_max, points)
    k = np.array([beb.maxwellian_integral(float(x))['rate_within_energy_range_m3_s'] for x in t])
    lt, lk = np.log(t), np.log(np.maximum(k, 1e-300))

    def rate(temperature):
        temp = np.asarray(temperature, float)
        out = np.exp(np.interp(np.log(np.clip(temp, t_min, t_max)), lt, lk))
        return np.where(temp < t_min, 0., out)
    return rate


def apply_ionization(mesh, conserved, entropy, magnetic, cold, neutral, dt, ion_mass, gamma, rate, cost_ev):
    rho = conserved[..., 0]
    momentum2 = conserved[..., 1]**2+conserved[..., 2]**2
    if conserved.shape[-1] == 6:
        momentum2 = momentum2+conserved[..., 4]**2/(.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None])
    kinetic_old = .5*momentum2/rho
    if entropy is not None:
        thermal = np.where(cold, entropy*rho**(gamma-1)/(gamma-1), conserved[..., 3]-kinetic_old-magnetic)
    else:
        thermal = conserved[..., 3]-kinetic_old-magnetic
    thermal = np.maximum(thermal, 0.)
    n_i = rho/ion_mass
    temperature = (gamma-1)*thermal/(2*n_i*QE)
    dn = neutral*-np.expm1(-n_i*rate(temperature)*dt)
    if cost_ev > 0:
        dn = np.minimum(dn, .5*thermal/(cost_ev*QE))
    add = dn*ion_mass
    new_rho = rho+add
    kinetic_new = .5*momentum2/new_rho
    cost = dn*cost_ev*QE
    result = conserved.copy()
    result[..., 0] = new_rho
    result[..., 3] -= cost
    new_entropy = entropy
    if entropy is not None:
        heat = thermal+(kinetic_old-kinetic_new)-cost
        new_entropy = np.where(dn > 0, (gamma-1)*np.maximum(heat, 0)*new_rho**(1-gamma), entropy)
    vol = mesh.volumes
    return result, new_entropy, neutral-dn, dict(mass_kg=float(np.sum(add*vol)), ions=float(np.sum(dn*vol)),
                                                 cost_j=float(np.sum(cost*vol)),
                                                 loading_heat_j=float(np.sum((kinetic_old-kinetic_new)*vol)))


class NeutralReservoirMHD(ConservativeVacuumMHD):
    """ConservativeVacuumMHD with a static neutral reservoir; attach with attach_reservoir after construction."""
    def attach_reservoir(self, neutral_m3, ion_mass, cost_ev, rate=None):
        self.neutral = np.asarray(neutral_m3, float).copy()
        self.neutral_initial = float(np.sum(self.neutral*self.mesh.volumes))
        self.reservoir_ion_mass, self.ionization_cost_ev = float(ion_mass), float(cost_ev)
        self.ionization_rate = beb_rate() if rate is None else rate
        self.ionized_ions = self.ionized_mass_kg = self.ionization_cost_j = self.loading_heat_j = 0.

    def _attempt(self, dt):
        result = super()._attempt(dt)
        if getattr(self, 'neutral', None) is None or not np.any(self.neutral > 0):
            return result
        state, p = self.state, self.parameters
        flux = self.mapping.fluid_flux(self.total_flux(state, self.time))
        cold = (self._cold(state.conserved, flux, state.entropy) if p.dual_energy
                else np.zeros(state.conserved.shape[:2], dtype=bool))
        conserved, entropy, self.neutral, booked = apply_ionization(
            self.mesh, state.conserved, state.entropy, self._cell_magnetic(state, flux), cold, self.neutral, dt,
            self.reservoir_ion_mass, p.gamma, self.ionization_rate, self.ionization_cost_ev)
        new = ConservativeState(conserved, state.induced_psi, entropy)
        self._primitive(new, self.time)
        self.state = new
        # The solver's mass and energy residuals subtract floor_mass/floor_energy as declared closure sources.
        self.floor_mass += booked['mass_kg']
        self.floor_energy -= booked['cost_j']
        self.ionized_ions += booked['ions']
        self.ionized_mass_kg += booked['mass_kg']
        self.ionization_cost_j += booked['cost_j']
        self.loading_heat_j += booked['loading_heat_j']
        return result

    def reservoir_diagnostics(self):
        left = float(np.sum(self.neutral*self.mesh.volumes))
        return dict(neutral_initial=self.neutral_initial, neutral_left=left,
                    neutral_left_fraction=left/self.neutral_initial if self.neutral_initial else 0.,
                    ionized_ions=self.ionized_ions, ionized_mass_kg=self.ionized_mass_kg,
                    ionization_cost_j=self.ionization_cost_j, loading_heat_j=self.loading_heat_j,
                    ionization_cost_ev=self.ionization_cost_ev)
