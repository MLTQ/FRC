"""Bounded compressible Hall MHD with compatible energy and implicit induction."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E
from .mhd import MHDParameters
from .toroidal_mhd import ToroidalMHD
from .toroidal_state import primitive, radius_squared
from .toroidal_flux import fast_speed
from .toroidal_reconstruction import ToroidalReconstruction
from .hall_compatible import CompatibleHallCurl
from .hall_midpoint import HallMidpoint
from .compatible_poynting import energy_rate


class HallMHD(ToroidalMHD):
    """Massless cold electrons, adiabatic ions; no vacuum or formation boundary."""
    def __init__(self, mesh, initial, ion_mass, parameters=None, *, hall=True, quadrature_order=5,
                 nonlinear_tolerance=1e-10, max_iterations=16):
        if not np.isfinite(ion_mass) or ion_mass <= 0 or type(hall) is not bool:
            raise ValueError('Positive finite ion mass and explicit boolean Hall choice required')
        if initial.magnetic_quadrature != 'compatible':
            raise ValueError('HallMHD requires explicitly constructed compatible magnetic energy; no energy conversion applied')
        super().__init__(mesh, initial, parameters or MHDParameters())
        self.ion_mass, self.hall_enabled = ion_mass, hall
        self.ideal_parameters = replace(self.parameters, resistivity=0.)
        self.geometry = CompatibleHallCurl(mesh, quadrature_order)
        self.induction = HallMidpoint(self.geometry, tolerance=nonlinear_tolerance, max_iterations=max_iterations)
        self.hall_iterations = 0
        self.max_hall_residual = self.max_linear_residual = 0.

    def timestep_limit(self):
        q = primitive(self.mesh, self.state, self.parameters.gamma)
        r = np.sqrt(radius_squared(self.mesh))
        ar = abs(q[..., 1])+fast_speed(q, r, 0, self.parameters.gamma)
        az = abs(q[..., 2])+fast_speed(q, r, 1, self.parameters.gamma)
        return self.parameters.cfl/np.max(ar/self.mesh.dr+az/self.mesh.dz)

    def density_samples(self, state):
        q = primitive(self.mesh, state, self.parameters.gamma)
        reconstruction = ToroidalReconstruction(self.mesh, q, order=self.parameters.spatial_order,
                                                 positivity=self.parameters.positivity_reconstruction)
        density = np.empty_like(self.geometry.r)
        for i in range(self.geometry.order):
            radius = self.geometry.r[:, 0, i, 0]
            for j in range(self.geometry.order):
                offset = (self.geometry.z[0, 0, 0, j]-self.mesh.z[0])/self.mesh.dz
                density[:, :, i, j] = reconstruction.raw(radius, offset)[..., 0]/self.ion_mass
        if np.any(density <= 0) or not np.isfinite(density).all():
            raise ValueError('Nonpositive carrier trace; no floor applied')
        return density

    def rhs(self, state):
        ideal = ToroidalMHD(self.mesh, state, self.ideal_parameters)
        derivative, _ = ideal.rhs(state)
        geometry = self.geometry
        current = geometry.current(state.psi, state.conserved[..., 5])
        electric = self.parameters.resistivity*current
        if self.hall_enabled:
            density = self.density_samples(state)
            magnetic = geometry.magnetic_samples(state.psi, state.conserved[..., 5])
            electric += np.cross(current, magnetic)/(density[..., None]*E)
        coefficients = geometry.project_electric(electric)
        dp, dg = geometry.unpack(geometry.faraday @ coefficients)
        derivative.psi += dp
        derivative.conserved[..., 5] += dg
        derivative.conserved[..., 3] += energy_rate(geometry, state.psi, state.conserved[..., 5], coefficients)
        heat = float(self.parameters.resistivity*np.sum(geometry.weight*np.sum(current*current, axis=-1)))
        return derivative, heat

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f'Time step must be positive and <= fluid limit {limit:.6g} s')
        first = ToroidalMHD(self.mesh, self.state, self.ideal_parameters)
        first.step(.5*dt)
        density = self.density_samples(first.state) if self.hall_enabled else None
        result = self.induction.step(first.state.psi, first.state.conserved[..., 5], density, dt,
                                     resistivity=self.parameters.resistivity)
        middle = first.state.copy()
        middle.psi = result.psi
        middle.conserved[..., 5] = result.g
        middle.conserved[..., 3] += result.total_energy_change
        primitive(self.mesh, middle, self.parameters.gamma)
        second = ToroidalMHD(self.mesh, middle, self.ideal_parameters)
        second.step(.5*dt)
        # Commit only after every magnetic, energy and fluid trial has passed.
        self.state = second.state
        self.time += dt
        self.steps += 1
        self.physical_ohmic_heat += result.physical_heat_j
        self.hall_iterations += result.iterations
        self.max_hall_residual = max(self.max_hall_residual, result.nonlinear_residual)
        self.max_linear_residual = max(self.max_linear_residual, result.linear_residual)
        for key in self.limiter_totals:
            if key == 'minimum_pressure_blend_weight':
                self.limiter_totals[key] = min(self.limiter_totals[key], first.limiter_totals[key], second.limiter_totals[key])
            else:
                self.limiter_totals[key] += first.limiter_totals[key]+second.limiter_totals[key]
        return dt

    def diagnostics(self):
        out = super().diagnostics()
        out.update(hall_enabled=self.hall_enabled, magnetic_energy_quadrature='compatible',
                   hall_nonlinear_iterations=self.hall_iterations, max_hall_nonlinear_residual=self.max_hall_residual,
                   max_hall_linear_residual=self.max_linear_residual, hall_matrix_factorizations=self.induction.factorizations)
        out['reconstruction_limiter']['accepted_cell_stages'] = 4*self.steps*self.mesh.nr*self.mesh.nz
        return out
