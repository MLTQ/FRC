"""Declared flux-limited parallel electron-conduction loss on open field lines; every joule is counted."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import e as QE, m_e as ME


@dataclass
class LossResult:
    conserved: np.ndarray
    entropy: np.ndarray | None
    energy_j: float
    cells: int


def open_line_loss(mesh, conserved, entropy, magnetic, cold, psi_cells, dt, ion_mass, gamma, flux_limit, wall_temperature_ev):
    """Remove thermal energy from cells on open field lines (psi >= 0, i.e. outside the separatrix).

    Electron conduction along open lines to the end walls is limited by free streaming,
    q <= f n T_e v_te with v_te = sqrt(e T_e/m_e). Over a connection length L (to the nearer
    tube end), the thermal energy (above the wall temperature) relaxes on
    tau = 1.5 L/(f v_te), integrated exactly as exp(-dt/tau) with tau frozen over the step.
    Te = Ti = p/(2 n e) as elsewhere; the ion share is lost with it (ends absorb both).
    Cold (entropy-pressure) cells update rho*A consistently.
    """
    rho = conserved[..., 0]
    kinetic = .5*(conserved[..., 1]**2+conserved[..., 2]**2)/rho
    if conserved.shape[-1] == 6:   # 2.5D: v_theta energy (the caller's `magnetic` includes B_theta energy)
        kinetic = kinetic+.5*conserved[..., 4]**2/(rho*.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None])
    thermal = np.where(cold, (entropy*rho**(gamma-1)/(gamma-1)) if entropy is not None else 0.,
                       conserved[..., 3]-kinetic-magnetic)
    n = rho/ion_mass
    te = np.maximum((gamma-1)*np.maximum(thermal, 0)/(2*n*QE), 0.)
    floor = n*2*QE*wall_temperature_ev/(gamma-1)
    open_ = (psi_cells >= 0) & (thermal > floor)
    if not np.any(open_):
        return LossResult(conserved, entropy, 0., 0)
    length = np.broadcast_to(mesh.half_length-abs(mesh.z)[None, :], rho.shape)+.5*mesh.dz
    tau = 1.5*length/(flux_limit*np.sqrt(QE*np.maximum(te, 1e-6)/ME))
    removed = np.where(open_, (thermal-floor)*(1-np.exp(-dt/tau)), 0.)
    result = conserved.copy()
    result[..., 3] -= removed
    new_entropy = entropy
    if entropy is not None:
        new_entropy = np.where(open_, (gamma-1)*np.maximum(thermal-removed, 0)*rho**(1-gamma), entropy)
    return LossResult(result, new_entropy, float(np.sum(removed*mesh.volumes)), int(open_.sum()))
