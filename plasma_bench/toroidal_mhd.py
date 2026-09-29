"""Conservative compressible r-z MHD with all three field/velocity components."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd import MHDParameters
from .mhd_induction import resistive_emf, resistive_energy_flux, ohmic_power
from .mhd_positivity import empty_statistics
from .toroidal_state import ToroidalState, primitive, integrated_quantities, radius_squared
from .toroidal_flux import fast_speed, rusanov
from .toroidal_reconstruction import ToroidalReconstruction
from .toroidal_induction import ideal_emf, toroidal_resistivity


class ToroidalMHD:
    """Bounded, fully filled reference cylinder; not the formation interface."""
    def __init__(self, mesh, initial, parameters=None):
        self.mesh, self.parameters = mesh, parameters or MHDParameters()
        primitive(mesh, initial, self.parameters.gamma)
        scale = max(float(np.max(abs(initial.psi))), 1e-30)
        if np.max(abs(initial.psi[0])) > 1e-13*scale or np.ptp(initial.psi[-1]) > 1e-13*scale:
            raise ValueError('Zero axis gauge and constant conducting radial-wall flux required')
        if not np.allclose(initial.psi[:, 0], initial.psi[:, -1], rtol=0., atol=1e-13*scale):
            raise ValueError('Periodic axial flux endpoints required')
        self.state = initial.copy()
        self.time, self.steps, self.physical_ohmic_heat = 0., 0, 0.
        self.limiter_totals = empty_statistics()
        self.initial_integrals = integrated_quantities(mesh, initial, self.parameters.gamma)

    def timestep_limit(self):
        q = primitive(self.mesh, self.state, self.parameters.gamma)
        r = np.sqrt(radius_squared(self.mesh))
        ar = abs(q[..., 1])+fast_speed(q, r, 0, self.parameters.gamma)
        az = abs(q[..., 2])+fast_speed(q, r, 1, self.parameters.gamma)
        wave = np.max(ar/self.mesh.dr+az/self.mesh.dz)
        diffusion = self.parameters.resistivity/MU0*(5/self.mesh.dr**2+2/self.mesh.dz**2)
        return self.parameters.cfl/(wave+diffusion)

    def _rhs(self, state):
        mesh, par = self.mesh, self.parameters
        q = primitive(mesh, state, par.gamma)
        reconstruction = ToroidalReconstruction(mesh, q, order=par.spatial_order, positivity=par.positivity_reconstruction)
        left, right, down, up = reconstruction.faces()
        br, bz = mesh.face_fields(state.psi)
        fr = rusanov(left, right, br, mesh.rf[:, None], 0, par.gamma)
        fz = rusanov(down, up, bz, np.sqrt(radius_squared(mesh)), 1, par.gamma)
        er = resistive_emf(mesh, state.psi, par.resistivity)
        sr, sz = resistive_energy_flux(mesh, q, er)
        dg, tr, tz, heat = toroidal_resistivity(mesh, q[..., 7], par.resistivity)
        fr[..., 3] += sr+tr
        fz[..., 3] += sz+tz
        du = -np.diff(mesh.rf[:, None, None]*fr, axis=0)/mesh.radial_measure[:, None, None]
        du -= np.diff(fz, axis=1)/mesh.dz
        du[..., 1] += reconstruction.radial_source()
        du[..., 5] += dg
        electric = ideal_emf(mesh, reconstruction, br, bz, par.gamma, par.spatial_order)+er
        return ToroidalState(du, -mesh.rf[:, None]*electric), heat+ohmic_power(mesh, er, par.resistivity), reconstruction.statistics

    def rhs(self, state):
        derivative, heat, _ = self._rhs(state)
        return derivative, heat

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f'Time step must be positive and <= {limit:.6g} s')
        old = self.state
        a, qa, sa = self._rhs(old)
        predictor = ToroidalState(old.conserved+dt*a.conserved, old.psi+dt*a.psi)
        primitive(self.mesh, predictor, self.parameters.gamma)
        b, qb, sb = self._rhs(predictor)
        new = ToroidalState(old.conserved+.5*dt*(a.conserved+b.conserved), old.psi+.5*dt*(a.psi+b.psi))
        primitive(self.mesh, new, self.parameters.gamma)
        self.state = new
        self.time += dt
        self.steps += 1
        self.physical_ohmic_heat += .5*dt*(qa+qb)
        for key in self.limiter_totals:
            if key == 'minimum_pressure_blend_weight':
                self.limiter_totals[key] = min(self.limiter_totals[key], sa[key], sb[key])
            else:
                self.limiter_totals[key] += sa[key]+sb[key]
        return dt

    def advance(self, end_time, *, max_steps=1_000_000, dt_max=np.inf):
        if not np.isfinite(end_time) or end_time < self.time or np.isnan(dt_max) or dt_max <= 0:
            raise ValueError('Finite nondecreasing end time and positive maximum step required')
        count = 0
        while self.time < end_time:
            remaining = end_time-self.time
            if remaining <= 4*np.finfo(float).eps*max(abs(end_time), 1e-300):
                self.time = end_time
                break
            if count >= max_steps:
                raise RuntimeError('Step budget exceeded; accepted partial state remains available')
            self.step(min(remaining, dt_max, self.timestep_limit()))
            count += 1
        return self.state

    def diagnostics(self):
        out = integrated_quantities(self.mesh, self.state, self.parameters.gamma)
        out.update(time_s=self.time, steps=self.steps, physical_ohmic_heat_j=self.physical_ohmic_heat)
        for name in ('mass_kg', 'total_energy_j', 'angular_momentum_kg_m2_s', 'toroidal_flux_wb'):
            out[name+'_change'] = out[name]-self.initial_integrals[name]
        out['reconstruction_limiter'] = {**self.limiter_totals, 'accepted_cell_stages': 2*self.steps*self.mesh.nr*self.mesh.nz}
        return out
