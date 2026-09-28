"""Developmental conservative r-z resistive MHD on a bounded fluid cylinder."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_state import MHDState, primitive, integrated_quantities
from .mhd_flux import ghost_cells, rusanov, fast_speed
from .mhd_induction import ideal_emf, resistive_emf, resistive_energy_flux, ohmic_power
from .mhd_reconstruction import Reconstruction


@dataclass(frozen=True)
class MHDParameters:
    gamma: float = 5/3
    resistivity: float = 0.0
    cfl: float = .25
    spatial_order: int = 2

    def __post_init__(self):
        if not np.isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must be finite and greater than one")
        if not np.isfinite(self.resistivity) or self.resistivity < 0:
            raise ValueError("Resistivity must be finite and nonnegative [ohm m]")
        if not np.isfinite(self.cfl) or not 0 < self.cfl <= .4:
            raise ValueError("Use 0 < CFL <= 0.4 for this explicit reference method")
        if isinstance(self.spatial_order, bool) or self.spatial_order not in (1, 2):
            raise ValueError("Spatial order must be 1 (reference) or 2 (limited reconstruction)")


class AxisymmetricMHD:
    """All cells contain fluid; no experimental plasma/vacuum interface yet."""
    def __init__(self, mesh, initial: MHDState, parameters=None):
        self.mesh = mesh
        self.parameters = parameters or MHDParameters()
        primitive(mesh, initial, self.parameters.gamma)
        scale = max(float(np.abs(initial.psi).max()), 1e-30)
        if np.max(np.abs(initial.psi[0])) > 1e-13*scale:
            raise ValueError("Axis flux gauge must be zero")
        if not np.allclose(initial.psi[:, 0], initial.psi[:, -1], rtol=0, atol=1e-13*scale):
            raise ValueError("Axial flux endpoints must be periodic")
        if np.ptp(initial.psi[-1]) > 1e-13*scale:
            raise ValueError("Slip conducting radial wall requires Br=0 (constant boundary psi)")
        self.state = initial.copy()
        self.time = 0.0
        self.steps = 0
        self.physical_ohmic_heat = 0.0
        self.initial_integrals = integrated_quantities(mesh, initial, self.parameters.gamma)

    def timestep_limit(self, state=None):
        q = primitive(self.mesh, self.state if state is None else state, self.parameters.gamma)
        ar = np.abs(q[..., 1])+fast_speed(q, 0, self.parameters.gamma)
        az = np.abs(q[..., 2])+fast_speed(q, 1, self.parameters.gamma)
        hyperbolic = float(np.max(ar/self.mesh.dr+az/self.mesh.dz))
        diffusion = 2*self.parameters.resistivity/MU0*(1/self.mesh.dr**2+1/self.mesh.dz**2)
        return self.parameters.cfl/(hyperbolic+diffusion)

    def rhs(self, state):
        mesh, par = self.mesh, self.parameters
        q = primitive(mesh, state, par.gamma)
        ghost = ghost_cells(q)
        brf, bzf = mesh.face_fields(state.psi)
        reconstruction = Reconstruction(mesh, q) if par.spatial_order == 2 else None
        if reconstruction is None:
            left, right, down, up = ghost[:-1, 1:-1], ghost[1:, 1:-1], ghost[1:-1, :-1], ghost[1:-1, 1:]
        else:
            left, right, down, up = reconstruction.faces()
        fr, _ = rusanov(left, right, brf, 0, par.gamma)
        fz, _ = rusanov(down, up, bzf, 1, par.gamma)
        er = resistive_emf(mesh, state.psi, par.resistivity)
        sr, sz = resistive_energy_flux(mesh, q, er)
        fr[..., 3] += sr
        fz[..., 3] += sz
        derivative = -np.diff(mesh.rf[:, None, None]*fr, axis=0)/mesh.radial_measure[:, None, None]
        derivative -= np.diff(fz, axis=1)/mesh.dz
        pt = q[..., 3]+(q[..., 4]**2+q[..., 5]**2)/(2*MU0)
        derivative[..., 1] += pt/mesh.r[:, None] if reconstruction is None else reconstruction.radial_source()
        emf = ideal_emf(mesh, ghost, brf, bzf, par.gamma, reconstruction)+er
        dpsi = -mesh.rf[:, None]*emf
        return MHDState(derivative, dpsi), ohmic_power(mesh, er, par.resistivity)

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f"Time step must be positive and <= {limit:.6g} s (wave + diffusion limit)")
        old = self.state
        a, qa = self.rhs(old)
        predictor = MHDState(old.conserved+dt*a.conserved, old.psi+dt*a.psi)
        primitive(self.mesh, predictor, self.parameters.gamma)
        b, qb = self.rhs(predictor)
        new = MHDState(old.conserved+.5*dt*(a.conserved+b.conserved), old.psi+.5*dt*(a.psi+b.psi))
        primitive(self.mesh, new, self.parameters.gamma)
        # Commit only after both stages validate; a failed step leaves state unchanged.
        self.state = new
        self.physical_ohmic_heat += .5*dt*(qa+qb)
        self.time += dt
        self.steps += 1
        return dt

    def advance(self, end_time, *, max_steps=1_000_000):
        if not np.isfinite(end_time) or end_time < self.time:
            raise ValueError("End time must be finite and not precede current state")
        count = 0
        while self.time < end_time:
            remaining = end_time-self.time
            if remaining <= 4*np.finfo(float).eps*max(abs(end_time), 1e-300):
                self.time = end_time
                break
            if count >= max_steps:
                raise RuntimeError("Step budget exceeded; state remains available for inspection")
            self.step(min(self.timestep_limit(), remaining))
            count += 1
        return self.state

    def diagnostics(self):
        out = integrated_quantities(self.mesh, self.state, self.parameters.gamma)
        out.update(time_s=self.time, steps=self.steps, physical_ohmic_heat_j=self.physical_ohmic_heat,
                   mass_relative_change=(out["mass_kg"]-self.initial_integrals["mass_kg"])/self.initial_integrals["mass_kg"],
                   total_energy_relative_change=(out["total_energy_j"]-self.initial_integrals["total_energy_j"])/self.initial_integrals["total_energy_j"])
        return out
