"""Dynamic vacuum mask for the 3D device hybrid (frc-ibl.11.9 prototype; see hybrid3d_vacuum_mask.md).

Gas cells holding fewer than `min_ions` macro-ions at both the start and the end of an ion step are treated as
exact vacuum for the field during that step. Their faces join the gap: B there is the potential field from
GapProjector (J = 0), and edges touching only masked or conductor cells carry no Ohm's-law E, current or Joule
heat. The true gas geometry (`types`) is unchanged for everything else: particle absorption, the gas-fraction
density normalization, and the electron energy closure.

Changing the mask re-projects the current field once at the start of the field step. With the boundary normal B
held, the potential field is the minimum-energy field in the newly masked region (Thomson's theorem), so this
removes any current that was flowing there and loses its energy. That loss is counted in
`mask_projection_response_j` and is not hidden in the load-energy residual.
"""
import time
import numpy as np
from .hybrid3d_device import DeviceHybrid, GapProjector, masks, GAS, GAP


class VacuumMaskHybrid(DeviceHybrid):
    def __init__(self, *args, min_ions=1, **kwargs):
        self.min_ions = int(min_ions)
        super().__init__(*args, **kwargs)
        if self.periodic_z or self.device.inner_radius > 0:
            raise NotImplementedError('Vacuum mask prototype: finite, non-annular tube only '
                                      '(a periodic or annular gap needs its harmonic fields rebuilt)')
        self.field_types = self.types.copy()
        self.mask_updates = 0
        self.mask_projection_response_j = 0.
        self.mask_projection_total_j = 0.
        self.mask_build_wall_s = 0.
        self.masked_cells = 0
        self._occupied_start = None
        self._apply_mask(self._occupancy())

    def _occupancy(self):
        """Cells holding at least min_ions macro-ions."""
        shape = tuple(int(s) for s in self.shape)
        idx = np.floor((self.x-self.lo)/self.d).astype(np.int64)
        for a in range(3):
            np.clip(idx[:, a], 0, shape[a]-1, out=idx[:, a])
        counts = np.bincount(np.ravel_multi_index(idx.T, shape), minlength=int(np.prod(shape)))
        return counts.reshape(shape) >= self.min_ions

    def _apply_mask(self, occupied):
        field = self.types.copy()
        field[(self.types == GAS) & ~occupied] = GAP
        if np.array_equal(field, self.field_types) and self.mask_updates:
            return
        start = time.time()
        response_before, total_before = self.response_magnetic_energy(), self.magnetic_energy()
        old = [self.face_active[a] | self.face_gap[a] for a in range(3)]
        self.field_types = field
        self.face_active, self.face_gap, self.edge_mask = masks(field)
        self.u8 = [m.astype(np.uint8) for m in self.edge_mask]
        self.f8 = [m.astype(np.uint8) for m in self.face_active]
        if self.source is not None:
            # Faces that become gap-conductor faces keep their present response (it is constant there:
            # conductor-surface edges carry the source's own inductive E).
            external = self.source(self.time)
            for a in range(3):
                fixed = old[a] & ~(self.face_active[a] | self.face_gap[a])
                self._fixed_response[a][fixed] = (self.b[a]-external[a])[fixed]
        self._gap_project = GapProjector(field, self.d)
        self.b = self.project(self.b, self.time)
        self.mask_projection_response_j += response_before-self.response_magnetic_energy()
        self.mask_projection_total_j += total_before-self.magnetic_energy()
        self.masked_cells = int(((self.types == GAS) & (field == GAP)).sum())
        self.mask_updates += 1
        self.mask_build_wall_s += time.time()-start

    def step(self):
        self._occupied_start = self._occupancy()
        return super().step()

    def field_substeps(self):
        # Called by DeviceHybrid.step after the push and deposit, before the field stages: the mask covers
        # cells empty at both ends of this ion step.
        occupied = self._occupancy()
        if self._occupied_start is not None:
            occupied |= self._occupied_start
        self._apply_mask(occupied)
        return super().field_substeps()

    def mask_diagnostics(self):
        return dict(masked_gas_cells=self.masked_cells, gas_cells=int((self.types == GAS).sum()),
                    mask_updates=self.mask_updates, mask_projection_response_j=self.mask_projection_response_j,
                    mask_projection_total_j=self.mask_projection_total_j, mask_build_wall_s=self.mask_build_wall_s)
