# hybrid3d_circuit.py

With gas geometry, the source basis is made discretely current-free in gas
by `VacuumCoilBasis`, including its compatible edge potential for reciprocity.

`edge_current` gives the known coil-field curl for source subtraction in
plasma Ohm's law. Total B still acts on particles; only plasma response
current enters Hall force and plasma Ohmic heating.

Full-shot 3D capacitor-bank boundary, using the existing switched-bank model
and actual coil geometry. No saved future 2D plasma field or current is used.
This is a new coupling path and requires a load/convergence comparison before
physical formation conclusions can be drawn.

- `CircuitDrive`: builds the same calibrated bias, main, re-boost and auxiliary
  banks as the source configuration and advances the bias phase to t=0.
  Charge calibration always spans at least 40 µs, including when a short
  verification replay ends before the main-current peak.
- `linkage`: Φ_s=∫J_plasma·A_s dV on the full Cartesian Yee grid. J is the
  curl of B minus the known coil source; conductor/image-current edges are
  excluded. Interface edge currents use their full dual volumes.
- `advance`: updates image current only in connected sections using
  dI_img=−L_connected⁻¹dΦ, then advances `CircuitStepper` and its continuous
  next-step source waveforms, as the 2D coupled formation path does.
- `snapshot`: actual branch currents, capacitor voltages, load work, linkage
  and the bank's own energy residual. This residual is not a total plasma
  energy balance.

Actual field-source work is W_source=−∫Phi·dI_source, measured directly with
ion-step trapezoidal quadrature of accepted linkage and the linear source
waveform. When a bank connects after plasma formation, Phi+L I_image has a
nonzero frozen reference. Bank `image_work_j` consequently differs from
field-source work by the reference-flux exchange, plus the measured source
waveform/circuit lag. `reference_flux_work_j`, `image_reference_flux_wb` and
`waveform_work_difference_j` retain these terms separately. This changes no
current, force, temperature or field. Comparing plasma load energy only with
bank image work falsely reported about11.4J loss in the end-plug case.

The start must be t=0; restoring a warm circuit state is unsupported. The
device source interface expects face B and edge ∂A/∂t at each field stage.
No azimuthal restriction is applied to plasma evolution or its circuit load.
