"""Declared minimum ionized density inside the chamber: a gap closure whose mass and energy are counted."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import e as QE


@dataclass
class FloorResult:
    conserved: np.ndarray
    entropy: np.ndarray | None
    mass_kg: float
    energy_j: float
    cells: int


def apply_floor(mesh, conserved, entropy, magnetic, cold, minimum_density, temperature_ev, ion_mass, gamma):
    """Raise rho to minimum_density by adding ions/electrons at rest at temperature_ev.

    Momentum is unchanged (new particles at rest), so the mass-loading kinetic-energy
    drop becomes heat inside the cell (total energy unchanged by it). The added
    particles' own thermal energy, 2 e T/(gamma-1) per ion, is the only energy source.
    Cold (entropy-pressure) cells update rho*A from their entropy thermal energy.
    """
    rho = conserved[..., 0]
    need = rho < minimum_density
    if not np.any(need):
        return FloorResult(conserved, entropy, 0., 0., 0)
    add = np.where(need, minimum_density-rho, 0.)
    new_rho = rho+add
    momentum2 = conserved[..., 1]**2+conserved[..., 2]**2
    if conserved.shape[-1] == 6:   # 2.5D: (rho v_theta)^2 = L^2/<r^2>; L and B_theta/r are unchanged
        momentum2 = momentum2+conserved[..., 4]**2/(.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None])
    kinetic_old, kinetic_new = .5*momentum2/rho, .5*momentum2/new_rho
    thermal_add = add/ion_mass*2*QE*temperature_ev/(gamma-1)
    result = conserved.copy()
    result[..., 0] = new_rho
    result[..., 3] += thermal_add
    new_entropy = entropy
    if entropy is not None:
        thermal = np.where(cold, entropy*rho**(gamma-1)/(gamma-1), conserved[..., 3]-kinetic_old-magnetic)
        thermal = thermal+(kinetic_old-kinetic_new)+thermal_add
        new_entropy = np.where(need, (gamma-1)*np.maximum(thermal, 0)*new_rho**(1-gamma), entropy)
    vol = mesh.volumes
    return FloorResult(result, new_entropy, float(np.sum(add*vol)), float(np.sum(thermal_add*vol)), int(need.sum()))
