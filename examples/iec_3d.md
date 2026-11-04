# iec_3d.py

## Purpose
First real full 3D electrostatic kinetic electron/ion experiment for frc-cac.2.
Runs and saves actual time, particles, self-consistent potential, electrode charge
and an uncorrected energy ledger. This is a collisionless development experiment,
not the historical air discharge or a validated POPS/fusion calculation.

## Model
Vlasov-Poisson: m dv/dt=q E, dx/dt=v; -epsilon0 Laplacian(phi)=rho. Actual electron
mass and configurable hydrogenic ion mass. Full open 15 cm cage, 2 mm wires,
explicitly assumed spherical 40 cm grounded chamber. No imposed central well,
neutrality closure, symmetry constraint, density floor or energy repair.
Initial spherical plasma seed is an input with independent random 3D species
positions and physical Maxwellians. Optional uniform wall electron source has
specified current, energy, normal direction; it is not gas-discharge physics.
Bias is V0(1-exp(-t/rise))*(1+depth*sin(2*pi*f*t)); rise=0 selects DC immediately.
No inferred resonance or desired frequency is embedded in the force.
`drive_mode=legacy` retains that waveform. New `dc`, `fixed` and `hops`
modes use `CageDrive` with explicit voltage response, sensor and warmup times.
HOPS (Harmonically Oscillating Plasma Sphere) observes old-state ion radius and
dilation speed to gate the positive pulse; no radius-frequency law is forced.

## Components and contracts
- `DEFAULTS`: declared exploratory parameters, not an optimum or hardware advice.
- `run`: new output directory only; immutable source copies/hashes at entry;
  conforming mesh, linear-FEM deposition/Poisson, kick-drift-kick full orbits.
- `voltage`: explicit DC/ramp/modulation drive. `atomic_json`: atomic progress.
- `ion_motion`/`ion_sensor` provide live-population RMS, anisotropy, centroid and
  coherent/residual energy diagnostics; absorption can bias that population.
- Saved control events include actual pulse intervals, never inferred resonance.
- `main`: `python -m examples.iec_3d OUT [--config JSON] [--mesh NPZ]`.
- Frame JSON stores actual positions/species/energies and a physical midplane
  potential cut. Up to 12,000 positions use an independent display RNG, which
  never changes source/loading randomness. Display count is independent of
  physical particle count; points are not persistent trajectory identities.
- `STOP` requests owned-run cancellation at a step boundary; partial data stay.
  A final actual-time frame is saved if stopping between scheduled frames.
- Each `state_NNNN.npz` retains full 3D nodal potential/charge and all particle
  positions/velocities/species/cells with exact time and physical macro weight.
  JSON is the lightweight viewer adapter; it is not the entire scientific state.
  Final particle/field NPZ files also include their actual physical timestamp.
- `saved_intervals` changes output cadence only. HOPS sources/controller and
  companion docs are included in provenance; a campaign executes frozen copies.
  A cached mesh must match the declared bulk/wire resolution. Legacy voltage
  arithmetic uses the original exact integer-step times.
- Optional `source_stop_ns` ends the ideal wall emitter at an exact integration
  boundary; omission preserves continuous emission and original arithmetic.
  Injected charge/energy stop accumulating, with no particle removal or reset.
  This permits source-free relaxation after a driven pulse. A single pulse can
  use `fixed` with zero warmup and a period longer than the entire run: its
  duty times period sets release, and no second on-command occurs in the run.

## Restart contract
Every saved state now writes an atomic complete `checkpoint.npz`, including
RNG/controller, full-step electric field, initial energy and all ledgers.
`run(... checkpoint=..., resume=True)` appends to the same experiment only with
identical configuration, geometry, physics hashes and complete matching history.
Original source/provenance stay; new restart source, provenance and prior manifest
are archived under `restarts/`. Existing scientific states cannot be overwritten.
Exact integer-step times and the original kick/drift/solve/kick loop are retained.
`stop_step` bounds an actual restart comparison and reports checkpointed rather
than complete when it ends before the configured physical endpoint. CLI flags
are `--checkpoint`, `--resume`, `--stop-step`. A legacy particle frame alone is
insufficient; its reconstructed state requires an exact original-interval replay.

## Budgets and limitations
Field energy=phi^T K phi/2; conductor charge=(K phi-q) on Dirichlet nodes.
Source work integrates V*dQ_supply, where dQ_supply=dQ_electrode-q_collected.
Grounded-wall emission contributes real injected kinetic energy; cage/wall hits
remove actual particle charge/kinetic energy at the first mesh-face crossing.
Emission occurs before the full step's first kick/drift, with no unmatched
half-kick at birth. Timestep counts are formed in ns/ps units to avoid an extra
step from floating-point unit conversion; the exact dt is saved in resolution.
Report K+field+absorbed-injected-sourcework minus its initial value, without
rescaling. Midpoint electrode work and boundary-hit leapfrog losses have finite
dt error. Tetra-constant E, planar steel, macroparticle self-force/noise and
initial/evolving Debye resolution need convergence. Numerical conservation does
not validate confinement. No neutral/Coulomb collisions, secondary emission,
ionization, radiation, reactions or fusion gain. A short electron-timescale run
cannot establish ion resonance or an optimum. Gmsh is an optional dependency.
