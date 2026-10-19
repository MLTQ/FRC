"""Conservative fluid state and an ideal-gas pressure recovery with no floors."""
from dataclasses import dataclass
import numpy as np
from .magnetic_quadrature import selected_energy


class NonphysicalState(ValueError):
    """Numerics or input produced an invalid state; never silently repaired."""


@dataclass
class MHDState:
    conserved: np.ndarray  # rho, rho vr, rho vz, total energy density
    psi: np.ndarray

    def copy(self):
        return MHDState(self.conserved.copy(), self.psi.copy())


def primitive(mesh, state, gamma, *, energy_quadrature="cell", pressure_override=None):
    """Primitive variables; `pressure_override=(mask, pressure)` substitutes declared cell pressures
    (dual-energy cold cells) before the positivity check."""
    u = state.conserved
    if u.shape != (mesh.nr, mesh.nz, 4) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(state.psi)):
        raise NonphysicalState("Conserved state/flux must be finite with the declared mesh shapes")
    rho = u[..., 0]
    if np.any(rho <= 0):
        raise NonphysicalState(f"Nonpositive density: minimum {rho.min():.6g} kg/m^3")
    vr, vz = u[..., 1]/rho, u[..., 2]/rho
    br, bz = mesh.cell_fields(state.psi)
    magnetic = selected_energy(mesh, state.psi, energy_quadrature)
    pressure = (gamma-1)*(u[..., 3]-.5*rho*(vr*vr+vz*vz)-magnetic)
    if pressure_override is not None:
        pressure = np.where(pressure_override[0], pressure_override[1], pressure)
    if not np.all(np.isfinite(pressure)) or np.any(pressure <= 0):
        raise NonphysicalState(f"Nonpositive pressure: minimum {pressure.min():.6g} Pa; no pressure floor applied")
    return np.stack((rho, vr, vz, pressure, br, bz), axis=-1)


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
