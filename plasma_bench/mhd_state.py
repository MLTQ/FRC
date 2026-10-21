"""Conservative fluid state and an ideal-gas pressure recovery with no floors."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .magnetic_quadrature import selected_energy


class NonphysicalState(ValueError):
    """Numerics or input produced an invalid state; never silently repaired."""


@dataclass
class MHDState:
    conserved: np.ndarray  # rho, rho vr, rho vz, total energy density [, L = r rho v_theta, g = B_theta/r]
    psi: np.ndarray

    def copy(self):
        return MHDState(self.conserved.copy(), self.psi.copy())


def radius_squared(mesh):
    """<r^2> over each annular cell with weight r dr: (r_lo^2+r_hi^2)/2, shape (nr, 1)."""
    return .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]


def toroidal_energies(mesh, conserved):
    """Cell kinetic and magnetic energy densities of v_theta and B_theta from L = rho <r^2> Omega and
    g = B_theta/r (exact for Omega and g uniform over the cell): (L^2/(2 rho <r^2>), <r^2> g^2/(2 mu0)).
    None for a four-variable (poloidal-only) state."""
    if conserved.shape[-1] != 6:
        return None
    r2 = radius_squared(mesh)
    return .5*conserved[..., 4]**2/(conserved[..., 0]*r2), r2*conserved[..., 5]**2/(2*MU0)


def primitive(mesh, state, gamma, *, energy_quadrature="cell", pressure_override=None):
    """Primitive variables; `pressure_override=(mask, pressure)` substitutes declared cell pressures
    (dual-energy cold cells) before the positivity check.

    A six-variable (2.5D) state returns two more columns, Omega = v_theta/r and g = B_theta/r."""
    u = state.conserved
    toroidal = u.shape[-1] == 6
    if (u.shape != (mesh.nr, mesh.nz, 6 if toroidal else 4) or not np.all(np.isfinite(u))
            or not np.all(np.isfinite(state.psi))):
        raise NonphysicalState("Conserved state/flux must be finite with the declared mesh shapes")
    rho = u[..., 0]
    if np.any(rho <= 0):
        raise NonphysicalState(f"Nonpositive density: minimum {rho.min():.6g} kg/m^3")
    vr, vz = u[..., 1]/rho, u[..., 2]/rho
    br, bz = mesh.cell_fields(state.psi)
    magnetic = selected_energy(mesh, state.psi, energy_quadrature)
    pressure = (gamma-1)*(u[..., 3]-.5*rho*(vr*vr+vz*vz)-magnetic)
    if toroidal:
        kinetic_t, magnetic_t = toroidal_energies(mesh, u)
        pressure = pressure-(gamma-1)*(kinetic_t+magnetic_t)
    if pressure_override is not None:
        pressure = np.where(pressure_override[0], pressure_override[1], pressure)
    if not np.all(np.isfinite(pressure)) or np.any(pressure <= 0):
        raise NonphysicalState(f"Nonpositive pressure: minimum {pressure.min():.6g} Pa; no pressure floor applied")
    if toroidal:
        return np.stack((rho, vr, vz, pressure, br, bz, u[..., 4]/(rho*radius_squared(mesh)), u[..., 5]), axis=-1)
    return np.stack((rho, vr, vz, pressure, br, bz), axis=-1)


def with_toroidal(mesh, conserved, omega=0., toroidal_over_radius=0.):
    """Append L = rho <r^2> Omega and g = B_theta/r to a four-variable state, adding their energy to E."""
    shape = (mesh.nr, mesh.nz)
    omega, g = [np.broadcast_to(np.asarray(v, dtype=float), shape) for v in (omega, toroidal_over_radius)]
    r2 = radius_squared(mesh)
    out = np.concatenate((conserved, (conserved[..., 0]*r2*omega)[..., None], g[..., None].copy()), axis=-1)
    out[..., 3] += .5*conserved[..., 0]*r2*omega**2+r2*g**2/(2*MU0)
    return out


def from_primitive(mesh, psi, density, pressure, vr=0.0, vz=0.0, gamma=5/3, *, energy_quadrature="cell"):
    if not np.isfinite(gamma) or gamma <= 1:
        raise ValueError("Adiabatic index must be finite and greater than one")
    shape = (mesh.nr, mesh.nz)
    rho, p, ur, uz = [np.broadcast_to(np.asarray(v, dtype=float), shape) for v in (density, pressure, vr, vz)]
    energy = p/(gamma-1)+.5*rho*(ur*ur+uz*uz)+selected_energy(mesh, psi, energy_quadrature)
    state = MHDState(np.stack((rho, rho*ur, rho*uz, energy), axis=-1), np.asarray(psi, dtype=float).copy())
    primitive(mesh, state, gamma, energy_quadrature=energy_quadrature)
    return state


def integrated_quantities(mesh, state, gamma, *, energy_quadrature="cell"):
    q = primitive(mesh, state, gamma, energy_quadrature=energy_quadrature)
    vol = mesh.volumes
    return {"mass_kg": float(np.sum(q[..., 0]*vol)),
            "axial_momentum_kg_m_s": float(np.sum(state.conserved[..., 2]*vol)),
            "total_energy_j": float(np.sum(state.conserved[..., 3]*vol)),
            "internal_energy_j": float(np.sum(q[..., 3]/(gamma-1)*vol)),
            "kinetic_energy_j": float(np.sum(.5*q[..., 0]*(q[..., 1]**2+q[..., 2]**2)*vol)),
            "magnetic_energy_j": float(np.sum(selected_energy(mesh, state.psi, energy_quadrature)*vol)),
            "min_density_kg_m3": float(q[..., 0].min()), "min_pressure_pa": float(q[..., 3].min()),
            "max_velocity_m_s": float(np.hypot(q[..., 1], q[..., 2]).max()),
            "max_face_divergence_t_m": float(np.abs(mesh.divergence(state.psi)).max())}
