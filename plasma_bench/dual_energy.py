"""Opt-in dual-energy pressure for cold/low-beta cells; corrections are counted, not hidden."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .magnetic_quadrature import selected_energy


def entropy_density(rho, pressure, gamma):
    """rho*A with A = p/rho^gamma: an advected density in adiabatic flow."""
    return pressure*rho**(1-gamma)


def _speeds(mesh, rho, vr, vz, internal, psi, gamma):
    """Local fast-speed bound sqrt(gamma p/rho + B^2/(mu0 rho)); negative internal energy counts as zero."""
    br, bz = mesh.cell_fields(psi)
    return np.sqrt(gamma*(gamma-1)*np.maximum(internal, 0)/rho+(br*br+bz*bz)/(MU0*rho))


def entropy_rate(mesh, q, entropy, heating, gamma):
    """First-order Rusanov transport of rho*A plus the physical heating source (gamma-1) rho^(1-gamma) q.

    Faces use the larger adjacent |v_n|+c_f. Reflecting wall, ends and axis carry
    zero flux. Positive for the solver's CFL bound; not exactly the mass flux, so
    a uniform A is preserved only to the first/second-order density mismatch.
    """
    rho, vr, vz, p = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    fast = np.sqrt(gamma*p/rho+(q[..., 4]**2+q[..., 5]**2)/(MU0*rho))
    lr = np.maximum(abs(vr[:-1])+fast[:-1], abs(vr[1:])+fast[1:])
    lz = np.maximum(abs(vz[:, :-1])+fast[:, :-1], abs(vz[:, 1:])+fast[:, 1:])
    fr = np.zeros((mesh.nr+1, mesh.nz))
    fz = np.zeros((mesh.nr, mesh.nz+1))
    fr[1:-1] = .5*(entropy[:-1]*vr[:-1]+entropy[1:]*vr[1:])-.5*lr*(entropy[1:]-entropy[:-1])
    fz[:, 1:-1] = .5*(entropy[:, :-1]*vz[:, :-1]+entropy[:, 1:]*vz[:, 1:])-.5*lz*(entropy[:, 1:]-entropy[:, :-1])
    divergence = np.diff(mesh.rf[:, None]*fr, axis=0)/mesh.radial_measure[:, None]+np.diff(fz, axis=1)/mesh.dz
    return -divergence+(gamma-1)*rho**(1-gamma)*heating


def cold_cells(mesh, conserved, psi, gamma, threshold, energy_quadrature="cell", shock_threshold=None):
    """Cold = internal/total energy below threshold and NOT a compressive shock cell.

    Optional shock detector (shock_threshold not None): D = -div(v) * min(dr, dz) / c_f
    above it, with positive total-energy internal energy, keeps the total-energy
    pressure. Negative internal energy is always cold.
    """
    rho = conserved[..., 0]
    vr, vz = conserved[..., 1]/rho, conserved[..., 2]/rho
    kinetic = .5*rho*(vr*vr+vz*vz)
    magnetic = selected_energy(mesh, psi, energy_quadrature)
    internal = conserved[..., 3]-kinetic-magnetic
    if shock_threshold is None:
        return internal < threshold*conserved[..., 3], internal, kinetic, magnetic
    rvr = np.pad(mesh.r[:, None]*vr, ((1, 1), (0, 0)))
    rvr[0], rvr[-1] = -rvr[1], 0.  # axis mirror; wall normal velocity zero at face (approximate ghost)
    vzp = np.pad(vz, ((0, 0), (1, 1)))
    vzp[:, 0], vzp[:, -1] = -vzp[:, 1], -vzp[:, -2]
    div = (rvr[2:]-rvr[:-2])/(2*mesh.dr*mesh.r[:, None])+(vzp[:, 2:]-vzp[:, :-2])/(2*mesh.dz)
    shock = (-div*min(mesh.dr, mesh.dz)/_speeds(mesh, rho, vr, vz, internal, psi, gamma) > shock_threshold) & (internal > 0)
    return (internal < threshold*conserved[..., 3]) & ~shock, internal, kinetic, magnetic


@dataclass
class Selection:
    conserved: np.ndarray
    entropy: np.ndarray
    cold: np.ndarray
    energy_change_j: float
    absolute_change_j: float
    max_cold_ratio: float     # entropy thermal energy / total, over cold cells
    disagreements: int        # cold cells whose entropy thermal exceeds 10*threshold of total


def select(mesh, conserved, psi, entropy, gamma, threshold, energy_quadrature="cell", shock_threshold=None):
    """Cold cells take p = rho*A*rho^(gamma-1); `conserved` returns the reset total energy (used only in
    reset mode). Other cells resynchronize rho*A from total energy, so shock heating enters the entropy."""
    cold, internal, kinetic, magnetic = cold_cells(mesh, conserved, psi, gamma, threshold, energy_quadrature,
                                                   shock_threshold)
    rho, total = conserved[..., 0], conserved[..., 3]
    thermal = entropy*rho**(gamma-1)/(gamma-1)
    result = conserved.copy()
    result[..., 3] = np.where(cold, thermal+kinetic+magnetic, total)
    synced = np.where(cold, entropy, (gamma-1)*np.maximum(internal, 0)*rho**(1-gamma))
    change = (result[..., 3]-total)*mesh.volumes
    ratio = np.where(cold, thermal/total, 0.)
    return Selection(result, synced, cold, float(change.sum()), float(abs(change).sum()),
                     float(ratio.max(initial=0.)), int(np.sum(cold & (ratio > 10*threshold))))
