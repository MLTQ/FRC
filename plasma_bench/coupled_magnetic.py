"""Midpoint magnetic work, Ohmic heating and elliptic vacuum constraints."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import diags, hstack
from scipy.sparse.linalg import spsolve
from .coupled_hydro import gas_primitive
from .magnetic_dissipation import shock_resistance


@dataclass
class MagneticUpdate:
    gas: np.ndarray
    induced: np.ndarray
    source_work: float
    ohmic_heat: float
    lorentz_work: float
    exchange_residual: float
    iterations: int
    nonlinear_residual: float
    numerical_heat: float = 0.
    maximum_numerical_resistivity: float = 0.


class MagneticStepper:
    """Adjoint force/advection pairing; not an energy reset after an MHD step."""
    def __init__(self, mapping, resistivity, gamma=5/3, *, tolerance=1e-11, max_iterations=40, method="direct", shock_diffusion=0.):
        if not np.isfinite(resistivity) or resistivity < 0:
            raise ValueError("Resistivity must be finite and nonnegative [ohm m]")
        if not 0 < tolerance < 1e-3 or max_iterations < 2:
            raise ValueError("Use a positive tight nonlinear tolerance and at least two iterations")
        self.mapping, self.resistivity, self.gamma = mapping, resistivity, gamma
        self.tolerance, self.max_iterations = tolerance, max_iterations
        if method not in ("direct", "picard"):
            raise ValueError("Midpoint method must be 'direct' or 'picard'")
        self.method = method
        if not np.isfinite(shock_diffusion) or not 0 <= shock_diffusion <= 2:
            raise ValueError("Use a declared numerical shock-diffusion coefficient between zero and two")
        self.shock_diffusion = shock_diffusion
        self.resistance = np.divide(resistivity*mapping.node_radius, MU0*mapping.node_area,
                                    out=np.zeros(mapping.node_count), where=mapping.active)

    def advection(self, total_flux):
        m = self.mapping
        dr, dz = np.gradient(total_flux, m.field.dr, m.field.dz, edge_order=2)
        return hstack((diags(m.grid_vector(dr))@m.interpolate_r,
                       diags(m.grid_vector(dz))@m.interpolate_z), format="csr")

    def step(self, gas, induced, external_before, external_after, dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Magnetic step duration must be finite and positive")
        m, tau = self.mapping, dt/2
        q = gas_primitive(m.fluid, gas, self.gamma)
        uold = m.interior(induced)
        ext0, ext1 = m.grid_vector(external_before), m.grid_vector(external_after)
        ext_mid = .5*(np.asarray(external_before)+np.asarray(external_after))
        delta_external = ext1-ext0
        velocity = np.concatenate((q[..., 1].ravel(), q[..., 2].ravel()))
        mass = np.tile((q[..., 0]*m.fluid.volumes).ravel(), 2)
        inverse_mass = diags(1/mass, format="csr") if self.method == "direct" else None
        midpoint = uold.copy()
        velocity_mid = velocity.copy()
        scale = max(np.linalg.norm(uold), np.linalg.norm(delta_external), 1e-30)
        for iteration in range(1, self.max_iterations+1):
            c = self.advection(m.full(midpoint)+ext_mid)
            numerical_r, maximum_eta = shock_resistance(m, q, m.full(midpoint)+ext_mid, velocity_mid, self.shock_diffusion, self.gamma)
            resistance = self.resistance+numerical_r
            if self.method == "direct":
                mobility = diags(resistance, format="csr")+tau*(2*np.pi/MU0)*(c@inverse_mass@c.T)
                matrix = m.active_diagonal+tau*(mobility@m.k)+m.vacuum_rows
                rhs = m.active*(uold-.5*delta_external)-tau*(c@velocity)
                candidate = spsolve(matrix.tocsc(), rhs)
            else:
                ku = m.k@midpoint
                vmid = velocity+tau*(2*np.pi/MU0)*(c.T@ku)/mass
                candidate = m.complete_vacuum(uold-.5*delta_external-tau*(c@vmid)-tau*resistance*ku)
            if not np.all(np.isfinite(candidate)):
                raise RuntimeError("Magnetic/vacuum solve returned nonfinite flux")
            scale = max(scale, np.linalg.norm(candidate))
            change = np.linalg.norm(candidate-midpoint)/scale
            midpoint = candidate
            c = self.advection(m.full(midpoint)+ext_mid)
            current_integral = m.k@midpoint
            gradient_energy = (2*np.pi/MU0)*current_integral
            velocity_mid = velocity+tau*(c.T@gradient_energy)/mass
            numerical_r, maximum_eta = shock_resistance(m, q, m.full(midpoint)+ext_mid, velocity_mid, self.shock_diffusion, self.gamma)
            resistance = self.resistance+numerical_r
            residual = m.active*(2*(midpoint-uold)+delta_external+dt*(c@velocity_mid)+dt*resistance*current_integral)
            nonlinear_residual = float(np.linalg.norm(residual)/scale)
            if change <= self.tolerance and nonlinear_residual <= 10*self.tolerance:
                break
        else:
            raise RuntimeError(f"Magnetic midpoint did not converge in {self.max_iterations} iterations; residual={nonlinear_residual:.3g}")
        final_u = 2*midpoint-uold
        final_velocity = 2*velocity_mid-velocity
        physical_nodes = dt*(2*np.pi/MU0)*self.resistance*current_integral**2
        numerical_nodes = dt*(2*np.pi/MU0)*numerical_r*current_integral**2
        heat_nodes = physical_nodes+numerical_nodes
        heat_cells = np.asarray(m.heat_to_cells@heat_nodes).reshape(q.shape[:2])
        vr, vz = np.split(final_velocity, 2)
        vr, vz = vr.reshape(q.shape[:2]), vz.reshape(q.shape[:2])
        final_gas = gas.copy()
        final_gas[..., 1], final_gas[..., 2] = q[..., 0]*vr, q[..., 0]*vz
        kinetic_change = .5*q[..., 0]*((vr-q[..., 1])*(vr+q[..., 1])+(vz-q[..., 2])*(vz+q[..., 2]))
        final_gas[..., 3] += kinetic_change+heat_cells/m.fluid.volumes
        gas_primitive(m.fluid, final_gas, self.gamma)
        source_work = -float(gradient_energy[m.active]@delta_external[m.active])
        lorentz_work = float(np.sum(kinetic_change*m.fluid.volumes))
        heat, numerical_heat = float(physical_nodes.sum()), float(numerical_nodes.sum())
        field_change = float(gradient_energy@(final_u-uold))
        return MagneticUpdate(final_gas, m.full(final_u), source_work, heat, lorentz_work,
                              field_change+lorentz_work+heat+numerical_heat-source_work, iteration, nonlinear_residual,
                              numerical_heat, maximum_eta)
