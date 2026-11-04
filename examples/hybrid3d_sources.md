# hybrid3d_sources.py

With gas geometry, `VacuumCoilBasis` makes each source current-free there
while retaining boundary-normal flux and a compatible edge potential. It
changes only the independently known external field, not plasma modes.

`edge_current` supplies the discrete curl of the known coil field. Ohm's law
subtracts it: coil currents lie outside gas, and their coarse grid stencil
must never be accelerated or heated as if it were plasma current.

Maps the initial poloidal flux and physical coil geometry to Cartesian Yee
edges, and supplies prescribed coil-current waveforms for 3D continuations.

- `flux_edges`: Aθ=ψ/r, regularized with ψ/r² on axis. The discrete curl
  preserves div B. Source-grid edge values are held beyond the source domain,
  inside the outer conductor. It never rebuilds the evolving plasma field.
  Flux is normalized before the cubic fit, with explicit solver tolerances:
  per-ampere coil fluxes are below scipy's default absolute spline tolerance,
  which would otherwise produce an identically zero source interpolant.
- `BankWaveform`: precomputes unit-current edge A and face B. Currents are
  linearly interpolated from saved vacuum-equivalent bank histories;
  `edge_derivative` gives ∂A/∂t for boundary Faraday. Out-of-range times fail.

`hybrid3d_frc` supplies source coil fluxes and the exact restart time.
`DeviceHybrid` subtracts only this known external field for its current-free
vacuum response, and uses −∂A/∂t outside the gas when advancing harmonic flux.
The plasma field is fully 3D and self-consistent. The drive is prescribed:
the old 2D bank/plasma load is not recoupled to the new 3D response. This
limitation must be present in run metadata and displayed descriptions.
