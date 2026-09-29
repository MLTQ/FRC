"""Conservative axisymmetric state including angular momentum and toroidal flux."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_state import NonphysicalState
from .magnetic_quadrature import selected_energy


@dataclass
class ToroidalState:
    conserved: np.ndarray  # rho, mr, mz, E, L=r rho vtheta, g=Btheta/r
    psi: np.ndarray

    def copy(self):
        return ToroidalState(self.conserved.copy(), self.psi.copy())


def radius_squared(mesh):
    return .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]


def primitive(mesh, state, gamma):
    u = state.conserved
    if u.shape != (mesh.nr, mesh.nz, 6) or state.psi.shape != (mesh.nr+1, mesh.nz+1) or not np.isfinite(u).all() or not np.isfinite(state.psi).all():
        raise NonphysicalState('Finite toroidal state with six conserved variables and corner flux required')
    rho = u[..., 0]
    if np.any(rho <= 0):
        raise NonphysicalState('Nonpositive density in toroidal state; no floor applied')
    vr, vz = u[..., 1]/rho, u[..., 2]/rho
    omega, g = u[..., 4]/(rho*radius_squared(mesh)), u[..., 5]
    br, bz = mesh.cell_fields(state.psi)
    kinetic = .5*rho*(vr**2+vz**2+radius_squared(mesh)*omega**2)
    magnetic = selected_energy(mesh, state.psi, 'cell')+radius_squared(mesh)*g**2/(2*MU0)
    p = (gamma-1)*(u[..., 3]-kinetic-magnetic)
    if not np.isfinite(p).all() or np.any(p <= 0):
        raise NonphysicalState(f'Nonpositive toroidal pressure (minimum {p.min():.6g}); no energy repair applied')
    return np.stack((rho, vr, vz, p, br, bz, omega, g), axis=-1)


def from_primitive(mesh, psi, density, pressure, vr=0., vz=0., omega=0., toroidal_over_radius=0., gamma=5/3):
    if not np.isfinite(gamma) or gamma <= 1:
        raise ValueError('Finite gamma greater than one required')
    shape = (mesh.nr, mesh.nz)
    rho, p, vr, vz, omega, g = [np.broadcast_to(np.asarray(v, float), shape) for v in
                              (density, pressure, vr, vz, omega, toroidal_over_radius)]
    r2 = radius_squared(mesh)
    e = p/(gamma-1)+.5*rho*(vr**2+vz**2+r2*omega**2)+selected_energy(mesh, psi, 'cell')+r2*g**2/(2*MU0)
    state = ToroidalState(np.stack((rho, rho*vr, rho*vz, e, rho*r2*omega, g), axis=-1), np.asarray(psi, float).copy())
    primitive(mesh, state, gamma)
    return state


def integrated_quantities(mesh, state, gamma):
    q = primitive(mesh, state, gamma)
    r2, volume = radius_squared(mesh), mesh.volumes
    kinetic = .5*q[..., 0]*(q[..., 1]**2+q[..., 2]**2+r2*q[..., 6]**2)
    toroidal = r2*q[..., 7]**2/(2*MU0)
    magnetic = selected_energy(mesh, state.psi, 'cell')+toroidal
    return {'mass_kg': float(np.sum(volume*q[..., 0])),
            'axial_momentum_kg_m_s': float(np.sum(volume*state.conserved[..., 2])),
            'angular_momentum_kg_m2_s': float(np.sum(volume*state.conserved[..., 4])),
            'toroidal_flux_wb': float(np.sum(volume*q[..., 7])/(2*np.pi)),
            'total_energy_j': float(np.sum(volume*state.conserved[..., 3])),
            'internal_energy_j': float(np.sum(volume*q[..., 3]/(gamma-1))),
            'kinetic_energy_j': float(np.sum(volume*kinetic)),
            'magnetic_energy_j': float(np.sum(volume*magnetic)),
            'toroidal_magnetic_energy_j': float(np.sum(volume*toroidal)),
            'min_density_kg_m3': float(q[..., 0].min()), 'min_pressure_pa': float(q[..., 3].min()),
            'max_velocity_m_s': float(np.sqrt(q[..., 1]**2+q[..., 2]**2+r2*q[..., 6]**2).max()),
            'max_face_divergence_t_m': float(np.abs(mesh.divergence(state.psi)).max())}
