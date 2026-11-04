# hybrid3d_energy_electrons.py

Conservative electron internal-energy closure, introduced to correct measured
full-shot pressure-work errors and selected by the device's thermal default.
Its formation path remains under verification. `EnergyElectrons` shares Spitzer resistivity and
edge Joule heating with `ThermalElectrons`; it evolves U=p/(γ−1), not entropy
carriers.

The equation is ∂tU+div[(U+p)u_e]=u_e·grad(p)+ηJ², with γ=5/3 and
u_e=u_i−J/(ne). This is the conservative pressure form of Eq.15 in
[Stanier, Chacón & Chen (2018/2019)](https://arxiv.org/abs/1803.07158).
Our explicit Yee/CIC implementation does not inherit that paper's implicit
midpoint conservation proof.

- `ion_pressure_work`: deposits the actual local Boris pressure-force power
  and work velocity, using old positions, mean pre/post velocities and the
  same CIC weights. Includes subsequently absorbed ions. This supplies a
  physical work reference, not a correction based on total energy error.
- `hall_pressure_work`: J_node·E_pressure, per gas dual volume. The caller
  must supply the same RK-stage weights as the field equation.
- `advance`: conservative upwind internal-energy flux and pressure compression
  on partial gas dual volumes, plus Joule heat. Compression is the negative
  adjoint of the central pressure gradient: −p div(f u), with gas fraction f.
  Its global work equals f u·grad(p), while it cannot extract heat from a
  zero-energy node. The link areas and cut-wall areas
  come from the same staircase gas cells. Outgoing wall enthalpy is counted;
  no incoming vacuum heat is supplied. The boundary flux includes enthalpy.
  Upwind transport is first order and
  diffusive. A materially negative energy throws instead of clipping it.
- `finish_step`: removes and counts energy on nodes with no actual PIC ions;
  no electron pressure or heat capacity remains in empty material.
- `transport_fluxes`: an optional face-flux/wall-rate provider. The CPU
  default returns None and retains its original calculations/order. The
  Metal subclass supplies float32-evaluated fluxes converted to float64.
  Both backends share the same conservative update, pressure work, numerical
  heat-removal prescription, positivity checks and energy counters.
  Providers also return the coefficient of outgoing upwind transport.
- Optional `fused_arrays=True` uses `hybrid3d_arrays.EnergyArrays` for
  float64 compression, CPU fluxes, pressure power, edge eta, Joule transfer
  and conservative RHS with fewer arrays. Stage work/heat results remain
  independently owned until their RK sum. The default remains NumPy; the
  verified fused path is selected explicitly. Metal keeps its
  own flux provider. The base closure still owns positivity checks, energy
  installation and every work/heat/boundary counter; physics is unchanged.

Electron transport has its own positivity timestep. With gas fraction f,
upwind outgoing rate A, wall rate b and central compression D(f u), set
lambda=(A+gamma*b+(gamma-1)*max(D(f u),0))/f on gas nodes. Advance the
electron operator in ceil(dt*max(lambda)/0.8) equal explicit substeps, with
velocity and stage heating held fixed. This keeps each node's own-U
coefficient nonnegative; all incoming transport and retained heat are
nonnegative. It includes wall enthalpy and expansion, not just one-axis CFL.
Pressure evolves during these split energy substeps; the associated explicit
field/ion pressure-centering discrepancy is measured, never fitted away.
The method remains first order/upwind and diffusive, with the actual energy
substep in the numerical diffusion estimate. Material negatives still throw.

This corrects an observed CPU 10 mm failure at 1.141 µs: an empty near-wall node
had positive U before its field substep but negative U after outgoing losses.
The old whistler/acoustic field bound did not guarantee transport positivity.
No energy/temperature clamp or symmetry constraint is introduced.
Independent periodic Gaussian translation at outer CFL 1.45–1.49 uses two
energy substeps and stays positive, with energy relative error ≤2.9e−16
and zero roundoff correction. L1 error decreases 0.825→0.591→0.387 over
32/64/128 cells; the narrow Gaussian is visibly diffused. The prescribed-flow
fixture is not formation. In the coupled 0.75 µs coarse replay, CPU and
Metal use at most two substeps, loss fraction ≤0.78221 and no roundoff
correction; grid-energy identity residual is about 1e−13 J, while the measured
pressure-centering residual is 0.00570 J. Full 10 mm shots remain underway.
- `sound_speed_squared`: uses the same density denominator as the force
  when setting the acoustic timestep; it does not change pressure or Te.

The existing Ohm-law density cutoff and its retained-heat fraction are still
declared approximations. Exact geometric exterior vacuum is preserved, with
no numerical atmosphere. Ionization, conduction and sheaths are absent.
Specifically, Joule heat retained is Q·(n/max(n,n_floor))^γ; the remainder
is explicitly removed. After deposition, U on nodes with n=0 is also removed
and counted. Neither is derived from a radiation or sheath model. Completed
v4 shots remove 527–780 J this way (43–60% of load work), and a combined
grid/time refinement at 0.75 µs changes those losses 20→106 J. These sinks
remain a central physical accuracy limitation, not a conservation success.
For constant 1D velocity, the upwind term has leading diffusion
D_num=|u| Δx(1−|u| Δt/Δx)/2; variable velocity, gas cuts and splitting add
errors beyond that simple estimate. See `docs/HYBRID_3D.md` for measured runs.
`absorbed_energy_j` now counts fluid wall enthalpy plus energy left on
emptied nodes; it is not a particle-carrier energy. Separate counters retain
both contributions. Stage energy residual and roundoff corrections are
measured without any global rescale of energy or temperature.

The caller must provide the actual Boris pressure work and RK Hall pressure
work/velocity. Explicit ion/field pressure time-centering error is measured
separately, without correcting the temperature. Independent contact/heating
comparisons and physical startup replays pass, and all six development shots
completed 40 µs. Formation/loss convergence remains unresolved; this class
must not be presented as a validated FRC confinement or lifetime model.
