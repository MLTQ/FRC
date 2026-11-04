"""Capacitor banks coupled to the full 3D plasma current by magnetic reciprocity."""
from dataclasses import replace
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.hybrid3d import curl_edges, curl_faces
from plasma_bench.coil_circuit import CircuitStepper
from examples.nucte_formation import coil_sections, circuit_drive
from examples.hybrid3d_sources import flux_edges
from plasma_bench.hybrid3d_coils import VacuumCoilBasis


class CircuitDrive:
    """I_coil=I_vac+I_image with image changes -L_connected^-1 dPhi_plasma.

    Phi_s = integral J_plasma . A_s dV. The full Cartesian current enters;
    neither it nor B is averaged over azimuth. Coils retain their specified
    axisymmetric physical geometry. This driver starts at the initial fill,
    not from a circuit state inferred at an arbitrary warm restart time.
    """
    def __init__(self, box, state, duration_s, gas_cells=None):
        if abs(state['time_s']) > 1e-12:
            raise ValueError('Coupled 3D circuit currently starts at t=0 only')
        project = VacuumCoilBasis(box, gas_cells) if gas_cells is not None else lambda a: a
        self.edges = np.array([project(flux_edges(box, state['grid'], psi)) for psi in state['unit_flux']])
        self.faces = np.array([curl_edges(*a, box.d) for a in self.edges])
        self.current_basis = np.array([np.asarray(curl_faces(*b, box.d))/MU0 for b in self.faces])
        self.coil = coil_sections(state['case'])
        nominal, self.calibration = circuit_drive(replace(state['case'], end_s=max(duration_s, 40e-6)), self.coil)
        self.stepper = CircuitStepper(self.coil.inductance, nominal.modules)
        self.stepper.run_to(0., 20e-9)
        self.previous_linkage = np.zeros(self.coil.sections)
        self.source_work_j = self.reference_flux_work_j = self.waveform_image_work_j = 0.

    def coefficients(self, time_s):
        return np.array([float(w.value(time_s)) for w in self.stepper.waveforms])

    def __call__(self, time_s):
        return np.einsum('s,saxyz->axyz', self.coefficients(time_s), self.faces)

    def edge_derivative(self, time_s):
        slopes = np.array([float(w.derivative(time_s)) for w in self.stepper.waveforms])
        return np.einsum('s,saxyz->axyz', slopes, self.edges)

    def edge_current(self, time_s):
        return np.einsum('s,saxyz->axyz', self.coefficients(time_s), self.current_basis)

    def linkage(self, model):
        response = np.asarray(model.b)-self(model.time)
        jp = np.asarray(curl_faces(*response, model.d))/MU0
        # Discrete edge currents already represent their full dual volumes,
        # including interface sheets. Exclude conductor/image current edges.
        jp *= np.asarray(model.edge_mask)
        return np.einsum('axyz,saxyz->s', jp, self.edges)*model.box.cell_volume

    def advance(self, model):
        phi = self.linkage(model)
        old_phi = self.previous_linkage.copy()
        old_image = self.stepper.image.copy()
        delta_current = self.coefficients(model.time)-self.coefficients(self.stepper.t)
        connected = self.stepper.connected()
        image = self.stepper.image.copy()
        if connected.any():
            image[connected] -= np.linalg.solve(self.coil.inductance[np.ix_(connected, connected)],
                                               (phi-self.previous_linkage)[connected])
        # A section connected after formation has a nonzero reference flux:
        # Phi + L I_image is constant while connected, not necessarily zero.
        # Measure actual source work directly; no state or energy is fitted.
        reference_old = old_phi+self.coil.inductance@old_image
        reference_new = phi+self.coil.inductance@image
        self.source_work_j -= float(.5*(old_phi+phi)@delta_current)
        self.reference_flux_work_j += float(.5*(reference_old+reference_new)@delta_current)
        self.waveform_image_work_j += float(.5*(old_image+image)@self.coil.inductance@delta_current)
        self.previous_linkage = phi
        self.stepper.advance(model.time, image)

    def snapshot(self):
        return dict(**self.stepper.snapshot(), plasma_linkage_wb=self.previous_linkage.tolist(),
                    source_work_j=self.source_work_j, reference_flux_work_j=self.reference_flux_work_j,
                    image_reference_flux_wb=(self.previous_linkage+self.coil.inductance@self.stepper.image).tolist(),
                    waveform_work_difference_j=self.waveform_image_work_j-self.stepper.snapshot()['image_work_j'],
                    circuit_energy_residual_j=self.stepper.energy_residual())
