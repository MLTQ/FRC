"""Conservative electron internal energy with local, compatible pressure work."""
import numpy as np
from scipy.constants import e as QE
from .hybrid3d import gradient_nodes
from .hybrid3d_electrons import ThermalElectrons


class EnergyElectrons(ThermalElectrons):
    """U=p/(gamma-1); U_t+div[(U+p)u_e]=u_e.grad(p)+eta J^2.

    Compression uses the negative adjoint of the pressure-force gradient;
    advection is conservative and upwind. Work velocities come from the
    actual Boris impulse and Faraday RK stages. No total-energy or
    temperature rescale is used. Under formation verification.
    """
    def __init__(self, model, gamma=5/3, resistivity='spitzer', constant_eta=0.,
                 coulomb_log=10., initial_temperature_ev=None, transport_cfl=.8,
                 fused_arrays=False):
        self.model, self.gamma = model, gamma
        self.resistivity_model, self.constant_eta = resistivity, constant_eta
        self.coulomb_log = coulomb_log
        self.fraction, self.d = model.node_gas_fraction, model.d
        self.live = self.fraction > 0
        self.floor, self.volume, self.q = model.p.n_floor, model.box.cell_volume, model.p.ion_charge
        te = model.p.te_ev if initial_temperature_ev is None else np.asarray(initial_temperature_ev)
        self.u = np.asarray(model.n*te*QE/(gamma-1)).copy()
        self.bulk_work = np.zeros_like(model.n)
        self.bulk_velocity = np.zeros((3,)+model.n.shape)
        gas = (model.types == 0).astype(float)
        self.link_area, self.boundary_plus, self.boundary_minus = [], [], []
        for a in range(3):
            p, q = [b for b in range(3) if b != a]
            def average(v):
                return .25*(v+np.roll(v, 1, p)+np.roll(v, 1, q)+np.roll(np.roll(v, 1, p), 1, q))
            self.link_area.append(average(gas))
            low = np.roll(gas, 1, a)
            self.boundary_plus.append(average(low*(1-gas))/self.d[a])
            self.boundary_minus.append(average(gas*(1-low))/self.d[a])
        self.ohmic_heat_j = self.cutoff_removed_heat_j = self.absorbed_energy_j = 0.
        self.enthalpy_outflow_j = self.vacuum_removed_energy_j = 0.
        self.max_transport_fraction = self.energy_balance_residual_j = 0.
        self.pressure_centering_residual_j = 0.
        self.roundoff_energy_j = 0.
        if not 0 < transport_cfl < 1:
            raise ValueError('Electron transport CFL must be between zero and one')
        self.transport_cfl = transport_cfl
        self.energy_substeps = 0
        self.max_energy_subcycles = 1
        self.max_energy_loss_fraction = 0.
        # Compatibility with the existing history schema; these operations
        # do not occur in this grid-energy formulation.
        self.remap_energy_j = self.heat_deposition_residual_j = 0.
        self.entropy_clip_count = 0
        self.arrays = None
        if fused_arrays:
            from .hybrid3d_arrays import EnergyArrays
            self.arrays = EnergyArrays(self)

    def pressure(self, n):
        return (self.gamma-1)*self.u

    def begin_step(self):
        self.bulk_work.fill(0.)

    def absorb(self, alive):
        # Electrons use the fluid enthalpy boundary flux, not particle IDs.
        pass

    def prepare_density(self, n):
        pass

    def ion_pressure_work(self, old_x, old_v, new_v, pressure_e):
        """Deposit the actual Boris pressure-work reference and work velocity.

        All pushed particles contribute, including those subsequently lost.
        CIC deposit is the adjoint of the force gather at the old positions.
        Compression uses its velocity; explicit pressure-centering error is
        measured against its power instead of correcting the temperature.
        """
        moment = self.model._moments(old_x, v=np.ascontiguousarray(.5*(old_v+new_v)))
        self.bulk_work = -self.q*np.sum(moment[1:]*pressure_e, axis=0)
        self.bulk_velocity = moment[1:]/np.maximum(self.model.n, self.floor)

    def hall_pressure_work(self, n, jn):
        if self.arrays is not None:
            return self.arrays.pressure_power(self.u, n, jn, self.gamma, self.q, self.floor)
        pressure_e = -gradient_nodes(self.pressure(n), self.d)/(self.q*np.maximum(n, self.floor))
        return np.sum(jn*pressure_e, axis=0)

    def edge_eta(self, n, masks):
        if self.arrays is not None:
            return self.arrays.edge_eta(self.eta(n), masks)
        return super().edge_eta(n, masks)

    def joule_density(self, currents, eta_edges):
        if self.arrays is not None:
            return self.arrays.heating(currents, eta_edges)
        return super().joule_density(currents, eta_edges)

    def transport_fluxes(self, u, velocity):
        """Optional (face fluxes, wall rate, outgoing rate); None keeps CPU."""
        return None

    def _transport_state(self, u, velocity):
        prepared = self.transport_fluxes(u, velocity)
        if prepared is not None:
            return prepared
        if self.arrays is not None:
            return self.arrays.fluxes(u, velocity)
        fields = []
        boundary_rate = np.zeros_like(u)
        outgoing = np.zeros_like(u)
        for a in range(3):
            speed = .5*(velocity[a]+np.roll(velocity[a], -1, a))
            fields.append(self.link_area[a]*speed*np.where(speed >= 0, u, np.roll(u, -1, a)))
            outgoing += (self.link_area[a]*np.maximum(speed, 0.)
                         +np.roll(self.link_area[a]*np.maximum(-speed, 0.), 1, a))/self.d[a]
            boundary_rate += (self.boundary_plus[a]*np.maximum(velocity[a], 0.)
                              +self.boundary_minus[a]*np.maximum(-velocity[a], 0.))
        return fields, boundary_rate, outgoing

    def advance(self, n, flux, jn, dt, heating, hall_work=None, hall_velocity=None):
        ne = np.maximum(n, self.floor)
        hall_velocity = -jn/(self.q*ne) if hall_velocity is None else hall_velocity
        velocity = self.bulk_velocity+hall_velocity
        if self.arrays is not None:
            compression = self.arrays.divergence_velocity(self.fraction, velocity)
        else:
            compression = np.zeros_like(n)
            for a in range(3):
                volume_velocity = self.fraction*velocity[a]
                compression += (np.roll(volume_velocity, -1, a)-np.roll(volume_velocity, 1, a))/(2*self.d[a])
        fields, boundary_rate, outgoing = self._transport_state(self.u, velocity)
        # Every coefficient removing a node's own U must fit within its gas
        # dual volume. The field whistler bound alone does not cover these
        # outgoing advective, enthalpy and expansion losses near cut walls.
        loss_rate = np.zeros_like(n)
        np.divide(outgoing+self.gamma*boundary_rate+(self.gamma-1)*np.maximum(compression, 0.),
                  self.fraction, out=loss_rate, where=self.live)
        loss_fraction = float(dt*np.max(loss_rate[self.live]))
        if not np.isfinite(loss_fraction):
            raise FloatingPointError('Nonfinite electron transport loss CFL')
        count = max(1, int(np.ceil(loss_fraction/self.transport_cfl)))
        if count > 1024:
            raise FloatingPointError('Electron transport needs >1024 substeps; reduce field/ion timestep')
        h = dt/count
        self.energy_substeps += count
        self.max_energy_subcycles = max(self.max_energy_subcycles, count)
        self.max_energy_loss_fraction = max(self.max_energy_loss_fraction, loss_fraction/count)
        self.max_transport_fraction = max(self.max_transport_fraction,
            float(np.max(np.abs(h*velocity[:, self.live]/self.d[:, None]))))
        retained = (np.maximum(n, 0.)/np.maximum(ne, 1.))**self.gamma
        heat = heating*retained
        hall_work = self.hall_pressure_work(n, jn) if hall_work is None else hall_work
        centered_work = (self.bulk_work+hall_work)*self.fraction
        for substep in range(count):
            before = self.energy(n)
            if substep:
                fields, _, _ = self._transport_state(self.u, velocity)
            if self.arrays is not None:
                full, work, outflow = self.arrays.rhs(self.u, self.fraction, fields,
                    boundary_rate, compression, heat, h, self.gamma)
            else:
                divergence = np.zeros_like(n)
                for a in range(3):
                    divergence += (fields[a]-np.roll(fields[a], 1, a))/self.d[a]
                outflow = self.gamma*self.u*boundary_rate
                work = -self.pressure(n)*compression
                full = self.u*self.fraction+h*(-divergence-outflow+work+heat)
            scale = max(float(np.max(np.abs(full))), 1.)
            if not np.isfinite(full[self.live]).all() or np.min(full[self.live]) < -1e-12*scale:
                index = np.unravel_index(np.argmin(np.where(self.live, full, np.inf)), full.shape)
                raise FloatingPointError(f'Nonpositive/nonfinite electron energy at t={self.model.time:g}, '
                    f'node={index}, loss_CFL={loss_fraction/count:g}; reduce timestep')
            self.pressure_centering_residual_j += float(h*(work-centered_work).sum()*self.volume)
            self.roundoff_energy_j += float(np.maximum(-full, 0.).sum()*self.volume)
            np.divide(np.maximum(full, 0.), self.fraction, out=self.u, where=self.live)
            self.u[~self.live] = 0.
            self.ohmic_heat_j += float(h*heating.sum()*self.volume)
            self.cutoff_removed_heat_j += float(h*(heating-heat).sum()*self.volume)
            escaped = float(h*outflow.sum()*self.volume)
            self.enthalpy_outflow_j += escaped
            self.absorbed_energy_j += escaped
            expected = float(h*(work+heat-outflow).sum()*self.volume)
            self.energy_balance_residual_j += self.energy(n)-before-expected

    def finish_step(self):
        empty = self.live & (self.model.n <= 0)
        removed = float(np.sum(self.u[empty]*self.fraction[empty])*self.volume)
        self.vacuum_removed_energy_j += removed
        self.absorbed_energy_j += removed
        self.u[empty] = 0.

    def energy(self, n):
        return float(np.sum(self.u*self.fraction)*self.volume)

    def sound_speed_squared(self, n):
        # The force denominator is the declared ne cutoff, so its acoustic
        # timestep uses that same denominator, not p/n in an empty cell.
        return self.gamma*self.pressure(n)/(self.model.p.ion_mass*np.maximum(n, self.floor))
