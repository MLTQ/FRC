"""Sharp plasma-vacuum edge for the 3D device hybrid (frc-ibl.11.9 prototype; design in hybrid3d_edge.md).

Plasma cells are gas cells holding at least `min_ions` macro-ions (at the start or the end of the ion step);
every other cell is exact vacuum for the field (VacuumMaskHybrid). In addition, Ohm's law at an interface node
compares the density of the adjacent plasma cells with the floor, not the density diluted over the node's
vacuum cells: with plasma-cell node fraction f_p and gas fraction f_g, the ideal and Hall terms use
ne = max(n, n_floor/s), s = min(f_g/f_p, max_scale), instead of max(n, n_floor). The cap keeps a node whose only
plasma neighbour holds a single, distant macro-ion (a CIC share near zero) from lowering the floor eightfold.
The ideal and Hall terms share ne, so the field-ion energy exchange stays exact. The electron-pressure term keeps
the closure's own density: max(n, n_floor) for EnergyElectrons, the interface density for EdgeEnergyElectrons.

With electron_edge (default), the energy closure gets the same edge (EdgeEnergyElectrons).
"""
import numpy as np
from scipy.constants import mu_0 as MU0
from .hybrid3d import gradient_nodes, curl_faces
from .hybrid3d_device import ohm_nodes_nb, ohm_edges_nb, curl_faces_nb, GAS
from .hybrid3d_energy_electrons import EnergyElectrons
from .hybrid3d_vacuum_mask import VacuumMaskHybrid


def node_fraction(cells):
    f = cells.astype(float)
    for a in range(3):
        f = .5*(f+np.roll(f, 1, axis=a))       # node (i,j,k) averages cells (i-1..i, j-1..j, k-1..k)
    return f


class EdgeEnergyElectrons(EnergyElectrons):
    """EnergyElectrons with the plasma-vacuum edge of PlasmaEdgeHybrid.

    1. Edge density. The pressure force on ions (in Ohm's law), the pressure work and the compression velocity
       n u/ne all use ne = max(n, n_floor/s), the model's interface density, instead of max(n, n_floor). They
       stay the same field and its adjoint, so the energy exchange remains exact; a diluted interface node no
       longer advects its electron energy n/n_floor too slowly behind a receding plasma.
    2. No internal-energy flux crosses a face that touches a vacuum node (no plasma cell around it).
       Quasi-neutral electrons cannot enter vacuum without ions; on a staircase edge, the first-order upwind
       flux of the tangential diamagnetic drift otherwise carries energy outward and never back.
    3. Energy on vacuum nodes is removed and counted at finish_step, as on emptied nodes.
    The update stays conservative; the conservation identity and positivity substeps are those of the base class.
    """
    def _plasma(self):
        nodes = getattr(self.model, '_plasma_nodes', None)
        return nodes if nodes is not None else self.model.n > 0

    def edge_density(self, n):
        s = getattr(self.model, '_scale', None)
        return np.maximum(n, self.floor) if s is None else np.maximum(n*s, self.floor)/s

    def ion_pressure_work(self, old_x, old_v, new_v, pressure_e):
        moment = self.model._moments(old_x, v=np.ascontiguousarray(.5*(old_v+new_v)))
        ne = self.edge_density(self.model.n)
        pressure_e = -gradient_nodes(self.pressure(self.model.n), self.d)/(self.q*ne)
        self.bulk_work = -self.q*np.sum(moment[1:]*pressure_e, axis=0)
        self.bulk_velocity = moment[1:]/ne

    def hall_pressure_work(self, n, jn):
        pressure_e = -gradient_nodes(self.pressure(n), self.d)/(self.q*self.edge_density(n))
        return np.sum(jn*pressure_e, axis=0)

    def _transport_state(self, u, velocity):
        fields, boundary_rate, outgoing = super()._transport_state(u, velocity)
        plasma = self._plasma()
        return [fields[a]*(plasma & np.roll(plasma, -1, a)) for a in range(3)], boundary_rate, outgoing

    def finish_step(self):
        empty = self.live & ((self.model.n <= 0) | ~self._plasma())
        removed = float(np.sum(self.u[empty]*self.fraction[empty])*self.volume)
        self.vacuum_removed_energy_j += removed
        self.absorbed_energy_j += removed
        self.u[empty] = 0.


class PlasmaEdgeHybrid(VacuumMaskHybrid):
    def __init__(self, *args, max_scale=2., electron_edge=True, **kwargs):
        self._scale, self.max_scale = None, float(max_scale)
        super().__init__(*args, **kwargs)
        if electron_edge and type(self.electrons) is EnergyElectrons and self.electrons.arrays is None:
            self.electrons.__class__ = EdgeEnergyElectrons
        self._zero_pressure = np.zeros(tuple(int(s) for s in self.shape))
        self.floor_bound_nodes = 0
        self.max_floor_ratio = 1.

    def _apply_mask(self, occupied):
        super()._apply_mask(occupied)
        fp = node_fraction(self.field_types == GAS)
        self._plasma_nodes = fp > 0
        self._scale = np.where(self._plasma_nodes,
                               np.minimum(self.node_gas_fraction/np.where(self._plasma_nodes, fp, 1.), self.max_scale), 1.)
        self._kernel_inv = self._scale*self._node_inv

    def electric(self, b, n, flux, time_s=None):
        if self._scale is None:                 # during DeviceHybrid.__init__, before the first mask
            return super().electric(b, n, flux, time_s)
        curl_faces_nb(b[0], b[1], b[2], *self.d, *self._j)
        if self.source is not None:
            t = self.time if time_s is None else time_s
            prescribed = (self.source.edge_current(t) if hasattr(self.source, 'edge_current') else
                          np.asarray(curl_faces(*self.source(t), self.d))/MU0)
            for a in range(3):
                self._j[a] -= prescribed[a]
        for a in range(3):
            self._j[a] *= self.edge_mask[a]
        s = self._scale
        # Scaling n, flux and J by s inside the kernel = gas-volume values with the floor lowered to n_floor/s.
        ohm_nodes_nb(b[0], b[1], b[2], *self._j, np.ascontiguousarray(n*s), np.ascontiguousarray(flux[0]*s),
                     np.ascontiguousarray(flux[1]*s), np.ascontiguousarray(flux[2]*s), self._zero_pressure,
                     self._kernel_inv, self.p.n_floor, self.p.ion_charge, *self.d, self._en, self._bn, self._jn)
        self._jn /= s                           # the electron closure and step() use gas-volume J
        pressure = self._pressure(n)
        if np.any(pressure):
            # the electron closure's own pressure density, so pressure work stays the adjoint of this force
            ne = (self.electrons.edge_density(n) if isinstance(self.electrons, EdgeEnergyElectrons)
                  else np.maximum(n, self.p.n_floor))
            self._en -= gradient_nodes(pressure, self.d)/(self.p.ion_charge*ne)
        ohm_edges_nb(self._en, *self._j, 0. if self.electrons is not None else self.p.eta,
                     self.p.eta_hyper, *self.d, *self.u8, *self._e)
        if self.electrons is not None:
            self._eta_edges = self.electrons.edge_eta(n, self.edge_mask)
            for a in range(3):
                self._e[a] += self._eta_edges[a]*self._j[a]
        return self._e

    def field_substeps(self):
        m = super().field_substeps()            # updates the mask, then the base bound with ne >= n_floor
        plasma = self._plasma_nodes
        if not plasma.any():
            return m
        ne = np.maximum(self.n[plasma]*self._scale[plasma], self.p.n_floor)/self._scale[plasma]
        ratio = self.p.n_floor/max(float(ne.min()), 1e-300)
        self.floor_bound_nodes = int(np.sum(self.n[plasma]*self._scale[plasma] < self.p.n_floor))
        self.max_floor_ratio = max(self.max_floor_ratio, ratio)
        # Whistler and Alfven rates grow at most like 1/ne: a lower interface ne scales the bound up.
        return max(m, int(np.ceil(m*ratio))) if ratio > 1 else m

    def mask_diagnostics(self):
        row = super().mask_diagnostics()
        row.update(floor_bound_plasma_nodes=self.floor_bound_nodes, max_interface_floor_ratio=self.max_floor_ratio)
        return row
