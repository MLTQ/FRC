"""Carrier-departure ledger and electron-edge transport on the production floor closure (frc-ibl.11.9).

CarrierLedgerHybrid is DeviceHybrid with its field, Ohm's law and density floor unchanged: no mask, no
projection, so no current is deleted. It adds:
- electron_edge (default on): the energy closure blocks internal-energy flux across any face touching a node with
  no ions at the start and the end of the ion step (EdgeEnergyElectrons). The model has no interface scale, so the
  closure's pressure force, pressure work and compression velocity are the base max(n, n_floor) exactly.
- a ledger that books, every RK stage with the solver's own weights, the Ohmic heat on edges by carrier class
  (from the macro-ions in the edge's four cells after the push): orphaned (none), sparse (fewer than
  `sparse_ions`, i.e. under one per cell on average) and resolved. It splits the cutoff-removed heat by node
  class (plasma n >= n_floor, thin, empty) and counts cells whose last macro-ion left.
The class sums reproduce the solver's ohmic and cutoff counters (checked in ledger_diagnostics).
"""
import numpy as np
from .hybrid3d_device import DeviceHybrid, GAS
from .hybrid3d_energy_electrons import EnergyElectrons
from .hybrid3d_edge import EdgeEnergyElectrons

EDGE_CLASSES = ('orphaned', 'sparse', 'resolved')
NODE_CLASSES = ('plasma', 'thin', 'empty')


class LedgerMixin:
    def ledger_init(self):
        self.ohmic_by_edge_class_j = dict.fromkeys(EDGE_CLASSES, 0.)
        self.cutoff_by_node_class_j = dict.fromkeys(NODE_CLASSES, 0.)
        self._stages = []

    def prepare_density(self, n):
        super().prepare_density(n)
        self.model._ledger_step(n)

    def joule_density(self, currents, eta_edges):
        heat = [eta_edges[a]*currents[a]**2 for a in range(3)]
        cls = self.model._edge_class
        self._stages.append([sum(float(np.sum(heat[a][cls[name][a]])) for a in range(3)) for name in EDGE_CLASSES])
        return super().joule_density(currents, eta_edges)

    def advance(self, n, flux, jn, dt, heating, hall_work=None, hall_velocity=None):
        q1, q2, q3 = self._stages[-3:]          # step() combines its three RK stages as (q1 + q2 + 4 q3)/6
        self._stages = []
        for i, name in enumerate(EDGE_CLASSES):
            self.ohmic_by_edge_class_j[name] += dt*(q1[i]+q2[i]+4*q3[i])/6*self.volume
        retained = (np.maximum(n, 0.)/np.maximum(np.maximum(n, self.floor), 1.))**self.gamma
        removed = heating*(1-retained)
        for name, m in zip(NODE_CLASSES, (n >= self.floor, (n > 0) & (n < self.floor), n <= 0)):
            self.cutoff_by_node_class_j[name] += dt*float(removed[m].sum())*self.volume
        return super().advance(n, flux, jn, dt, heating, hall_work=hall_work, hall_velocity=hall_velocity)


class LedgerElectrons(LedgerMixin, EnergyElectrons):
    pass


class EdgeLedgerElectrons(LedgerMixin, EdgeEnergyElectrons):
    pass


class CarrierLedgerHybrid(DeviceHybrid):
    def __init__(self, *args, electron_edge=True, sparse_ions=4, **kwargs):
        super().__init__(*args, **kwargs)
        e = self.electrons
        if type(e) is not EnergyElectrons or e.arrays is not None:
            raise ValueError('Carrier ledger needs the CPU EnergyElectrons closure without fused arrays')
        e.__class__ = EdgeLedgerElectrons if electron_edge else LedgerElectrons
        e.ledger_init()
        self.electron_edge, self.sparse_ions = bool(electron_edge), int(sparse_ions)
        self.gas = self.types == GAS
        self._counts = self._cell_counts()
        self._plasma_nodes = self.n > 0
        self.departed_cells = 0
        self._ledger_step(self.n)
        self.departed_cells = 0

    def _cell_counts(self):
        shape = tuple(int(s) for s in self.shape)
        idx = np.floor((self.x-self.lo)/self.d).astype(np.int64)
        for a in range(3):
            np.clip(idx[:, a], 0, shape[a]-1, out=idx[:, a])
        return np.bincount(np.ravel_multi_index(idx.T, shape), minlength=int(np.prod(shape))).reshape(shape)

    def _ledger_step(self, n_new):
        counts = self._cell_counts()
        # EdgeEnergyElectrons reads this: nodes with ions at the start or the end of the step
        self._plasma_nodes = (self.n > 0) | (n_new > 0)
        cls = {name: [] for name in EDGE_CLASSES}
        for a in range(3):
            p, q = [b for b in range(3) if b != a]
            c = counts+np.roll(counts, 1, p)+np.roll(counts, 1, q)+np.roll(np.roll(counts, 1, p), 1, q)
            cls['orphaned'].append(c == 0)
            cls['sparse'].append((c > 0) & (c < self.sparse_ions))
            cls['resolved'].append(c >= self.sparse_ions)
        self._edge_class = cls
        self.departed_cells += int(np.sum(self.gas & (self._counts > 0) & (counts == 0)))
        self._counts = counts

    def ledger_diagnostics(self):
        e = self.electrons
        area = [self.d[1]*self.d[2], self.d[0]*self.d[2], self.d[0]*self.d[1]]
        current = {name: float(sum(np.sum(np.abs(self._j[a])[self._edge_class[name][a]])*area[a] for a in range(3)))
                   for name in EDGE_CLASSES}
        return dict(ledger_electron_edge=self.electron_edge,
                    ohmic_by_edge_class_j=dict(e.ohmic_by_edge_class_j),
                    ohmic_class_sum_minus_counter_j=sum(e.ohmic_by_edge_class_j.values())-e.ohmic_heat_j,
                    cutoff_by_node_class_j=dict(e.cutoff_by_node_class_j),
                    cutoff_class_sum_minus_counter_j=sum(e.cutoff_by_node_class_j.values())-e.cutoff_removed_heat_j,
                    edge_current_magnitude_a=current, departed_cells=self.departed_cells,
                    occupied_gas_cells=int(np.sum(self.gas & (self._counts > 0))))
