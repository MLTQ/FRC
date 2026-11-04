"""Advected electron entropy and resistive heating for the kinetic-ion device model."""
import numpy as np
from scipy.constants import e as QE
from scipy.ndimage import distance_transform_edt, map_coordinates
from .hybrid3d import edges_to_nodes


class ThermalElectrons:
    """p=S n^gamma, DS/Dt=(gamma-1) eta J²/n^gamma, u_e=u_i-J/(ne).

    Bulk entropy is carried by the same macro-particles that supply density.
    Only relative Hall transport is remapped on the grid. Joule energy is
    distributed with the PIC weights and the exact discrete energy gradient.
    """
    def __init__(self, model, gamma=5/3, resistivity='spitzer', constant_eta=0., coulomb_log=10., initial_temperature_ev=None):
        self.gamma, self.resistivity_model = gamma, resistivity
        self.model = model
        self.constant_eta, self.coulomb_log = constant_eta, coulomb_log
        self.fraction = model.node_gas_fraction
        self.live = self.fraction > 0
        self.floor, self.volume = model.p.n_floor, model.box.cell_volume
        self.d, self.q, self.periodic_z = model.d, model.p.ion_charge, model.periodic_z
        self.periodic_axes = (not np.isfinite(model.device.wall_radius),)*2+(model.periodic_z,)
        self.index = np.indices(model.box.shape, dtype=float)
        self.nearest = distance_transform_edt(~self.live, return_distances=False, return_indices=True)
        initial_te = model.p.te_ev if initial_temperature_ev is None else np.asarray(initial_temperature_ev)
        self.entropy = initial_te*QE*np.maximum(model.n, 1.)**(1-gamma)
        self.entropy = self.entropy[tuple(self.nearest)]
        self.carried = self._sample(self.entropy)
        self.material_nearest = self.nearest
        self.entropy = self._grid_from_carriers()
        self.ohmic_heat_j = self.cutoff_removed_heat_j = self.max_transport_fraction = 0.
        self.remap_energy_j = self.heat_deposition_residual_j = self.absorbed_energy_j = 0.
        self.entropy_clip_count = 0

    def _sample(self, field):
        coordinates = (self.model.x.T-self.model.lo[:, None])/self.d[:, None]
        return map_coordinates(field, coordinates, order=1, mode='grid-wrap', prefilter=False)

    def _extend(self, field, n):
        positive = self.live & (n > 0)
        if self.material_nearest is None:
            return np.zeros_like(field)
        field = field.copy()
        field[~positive] = field[tuple(self.material_nearest)][~positive]
        return field

    def _grid_from_carriers(self):
        values = np.zeros_like(self.model.v)
        values[:, 0] = self.carried
        moments = self.model._moments(self.model.x, v=values)
        field = np.zeros_like(moments[0])
        np.divide(moments[1], moments[0], out=field, where=moments[0] > 0)
        return self._extend(field, moments[0])

    def begin_step(self):
        self.old_entropy = self.entropy.copy()
        self.old_carrier_energy = (self.model.w*self.carried*self._sample(np.maximum(self.model.n, 0.)**(self.gamma-1))
                                   /(self.gamma-1))

    def absorb(self, alive):
        self.absorbed_energy_j += float(self.old_carrier_energy[~alive].sum())
        self.carried = self.carried[alive]

    def prepare_density(self, n):
        positive = self.live & (n > 0)
        self.material_nearest = (distance_transform_edt(~positive, return_distances=False, return_indices=True)
                                 if positive.any() else None)
        self.bulk_end = self._grid_from_carriers()
        # Lie split: the PIC carriers supply bulk transport once per ion
        # step, then positive interpolation transports entropy relative to
        # them with the Hall velocity. Adding a fixed bulk-grid difference
        # repeatedly after Hall remaps could make entropy negative.
        self.entropy = self.bulk_end.copy()
        self.heat_energy = np.zeros_like(n)

    def pressure(self, n):
        return self.entropy*np.maximum(n, 0.)**self.gamma

    def temperature(self, n):
        return self.pressure(n)/(np.maximum(n, 1.)*QE)

    def eta(self, n):
        if self.resistivity_model == 'constant':
            eta = np.full_like(n, self.constant_eta)
        elif self.resistivity_model == 'spitzer':
            # NRL transverse Spitzer, SI; same lnLambda=10 and 1 eV
            # applicability cutoff as the formation solver. No T floor in p.
            eta = 1.03e-4*self.coulomb_log*np.maximum(self.temperature(n), 1.)**-1.5
        else:
            raise ValueError('Unknown electron resistivity model')
        return eta*self.live

    def edge_eta(self, n, masks):
        eta = self.eta(n)
        return [.5*(eta+np.roll(eta, -1, axis=a))*masks[a] for a in range(3)]

    def joule_density(self, currents, eta_edges):
        # Full dual-volume heating. The sum over nodes equals sum eta J²
        # over edges, including interface sheets, before a density cutoff.
        return np.sum(edges_to_nodes(*[eta_edges[a]*currents[a]**2 for a in range(3)]), axis=0)

    def advance(self, n, flux, jn, dt, heating):
        ne = np.maximum(n, self.floor)
        # The ion bulk motion is already represented by particle positions.
        relative_velocity = -jn/(self.q*np.maximum(ne, 1.))
        shift = dt*relative_velocity/self.d[:, None, None, None]
        self.max_transport_fraction = max(self.max_transport_fraction, float(np.max(np.abs(shift[:, self.live]))))
        departure = self.index-shift
        for a in range(3):
            if not self.periodic_axes[a]:
                departure[a] = np.clip(departure[a], 0., n.shape[a]-1.)
        entropy = map_coordinates(self.entropy, departure, order=1, mode='grid-wrap', prefilter=False)
        self.entropy = self._extend(entropy, n)
        self.heat_energy += dt*heating
        self.ohmic_heat_j += float(dt*np.sum(heating)*self.volume)
        if not np.isfinite(self.entropy[self.live]).all():
            raise FloatingPointError('Nonfinite electron entropy')

    def finish_step(self):
        n = self.model.n
        target = self.energy(n)
        delta = self.entropy-self.bulk_end
        decrease = np.zeros_like(delta)
        np.divide(np.maximum(-delta, 0.), self.bulk_end, out=decrease, where=self.bulk_end > 0)
        # Transfer a loss as a fraction of the existing carrier entropy,
        # and a gain as an absolute increment. Both are positive and their
        # weighted sums equal the requested nodal entropy changes. Direct
        # additive transfer could create negative carriers, whose clipping
        # would then manufacture heat.
        updated = self.carried*(1-self._sample(decrease))+self._sample(np.maximum(delta, 0.))
        self.entropy_clip_count += int(np.count_nonzero(updated < 0))
        self.carried = np.maximum(updated, 0.)
        self.entropy = self._grid_from_carriers()
        self.remap_energy_j += self.energy(n)-target
        retained = (np.maximum(n, 0.)/np.maximum(np.maximum(n, self.floor), 1.))**self.gamma
        heat = self.heat_energy*retained
        self.cutoff_removed_heat_j += float(np.sum(self.heat_energy-heat)*self.volume)
        per_ion = np.zeros_like(n)
        np.divide(heat, n*self.fraction, out=per_ion, where=(n > 0)&self.live)
        # Ue = sum_p w_p S_p gather[n^(gamma-1)]/(gamma-1).
        # Distribute local heat in proportion to the very same CIC weights;
        # this is an exact discrete energy increment, with no global rescale.
        coefficient = self._sample(np.maximum(n, 0.)**(self.gamma-1))
        increment = np.zeros_like(self.carried)
        np.divide((self.gamma-1)*self._sample(per_ion), coefficient, out=increment, where=coefficient > 0)
        before = self.energy(n)
        self.carried += increment
        self.entropy = self._grid_from_carriers()
        self.heat_deposition_residual_j += self.energy(n)-before-float(heat.sum()*self.volume)

    def energy(self, n):
        return float(np.sum(self.pressure(n)*self.fraction)*self.volume/(self.gamma-1))
