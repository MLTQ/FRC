# Reacting species with electric-field feedback

The full retained atomic moment state now drives an evolving electric field,
which in turn accelerates each charged species through its physical charge and
mass. Neutrals and exactly empty species receive no Lorentz force. The local
multifluid–Maxwell equations, fixed-frame transformation, energy/field-momentum
budgets and source splitting are in `plasma_bench/atomic_electromagnetic.md`.
The numerical source approach is based on
[Wang et al. (2020)](https://arxiv.org/abs/1909.04125).

This is a homogeneous physical source component. Displacement current and
electric energy are retained, unlike the moving solver's quasineutral constraint.
It is not installed in the moving pulse. Joining spatial species fluxes, Maxwell
curls, Gauss continuity and the exact moving vacuum interface remains required;
one cannot simply graft this source into the old fixed-carrier state.

## Complete midpoint histories

Run56187 exited0 in477.414s. All six analytic-mode histories and both full
700/1400-step heated argon histories through3.5us completed. Every accepted
species/E state and source are in `outputs/atomic_electromagnetic`, source
`5e55b33a3eab66fafc90f1558b62f2a65b913cff77c20f1cecfadea239e0f80e`.
The scientific PNG was inspected; PDF saved. No unit suite, floors, particle
normalization or energy repair was used.

Uniform electron/argon modes use physical masses, n=1e16/m³, Bz=20mT,
E0=100V/m and a fixed integration frame different from the laboratory frame.
Closed-form parallel and perpendicular mode references are specified in
`examples/atomic_electromagnetic.md`; they do not call another time integrator.

| Mode | Steps | Maximum E error [V/m] | Maximum J error [A/m²] |
|---|---:|---:|---:|
| Parallel | 100 | 29.757960 | 1.427425 |
| Parallel | 200 | 7.575867 | 0.361527 |
| Parallel | 400 | 1.904078 | 0.090606 |
| Perpendicular | 100 | 12.642758 | 0.819714 |
| Perpendicular | 200 | 3.201396 | 0.207769 |
| Perpendicular | 400 | 0.802260 | 0.052083 |

Errors fall by approximately four under timestep halving. Maximum energy drift
is7.11e-15eV/reference nucleus; total particle+field momentum drift<=1.93e-16
in the stated atomic momentum units. Thermal temperatures stay constant within
2.0e-15eV. These conservation results do not remove the measured phase error.

The first runner incorrectly included the undefined empty-neutral velocity in
one aggregate velocity-error metric. The original state/summary are retained;
`occupied_velocity_readout.json` recomputes the physical occupied-species values
from saved arrays. Parallel errors are890.916/225.644/56.551m/s; perpendicular
511.620/129.678/32.507m/s. The current driver excludes that empty slot. No field,
current, thermal or conservation result changes.

The complete heated argon cases use the prior joint excitation/ionization/
recombination/radiation/collision model and20eV/nucleus heat pulse, now with E
feedback and constant Bz=20mT. The finer final electron number/temperature are
0.683265159 and2.502599268eV. Maximum combined particle+electric+binding+photon
energy drift is2.66e-12eV/nucleus, charge error1.30e-13 and total scaled-momentum
drift2.06e-13. Photons remain produced energy, not proven escaped radiation.

**These ns-step midpoint current histories are not temporally resolved.**
Maximum omega_p*dt is2207.23/1103.62; cyclotron*dt is17.5882/8.7941. Time halving
changes E by up to3316.75V/m and J by27670.94A/m², although population L1 and Te
differences are only4.24815e-6 and0.000260378eV. The finer final current magnitude
is about19.70A/m², versus1.465A/m² coarse. Neither is a validated current tail.
This is concrete evidence that a small conserved-energy residual cannot certify
field/current dynamics.

## Completed exponential source histories

The same physical skew source can be integrated by a sparse matrix exponential,
with no midpoint phase dispersion inside a frozen-population substep. Optional
`method="exponential"` uses the actual physical frequencies and retains every
nonzero charged population, including trace states. Details and remaining
splitting limits are in `plasma_bench/species_field_exponential.md`.

Run77504 exited 0 after 566.191 s in `outputs/atomic_electromagnetic_exponential`, source
`8bae70d9b3335e0c0426c8978142fdc5bbfc77e931d46418ec4cf832a0052822`.
Its six complete analytic histories have maximum E error<=5.68e-13V/m,
J error<=3.50e-14A/m² and occupied velocity error<=2.18e-11m/s. Independent
energy drift<=7.11e-15eV/nucleus and total scaled-momentum drift<=2.07e-16.
The full700-step reacting history now completes3.5us: final current magnitude
0.008561895A/m², E=[-0.000447424,-2.40530e-7,-0.002872067]V/m, electron
number0.6832656324, Te2.502600599eV and energy drift2.52e-12eV/nucleus.
The full 1400-step trajectory also completes through 3.5 us. Its final current
magnitude is 0.003172620 A/m², electric field
[-0.000447383,-1.02326e-7,0.001229286] V/m, electron number 0.6832651966 and
Te=2.502599254 eV. Maximum energy-budget drift is 2.43e-12 eV/nucleus. The
scientific PNG was inspected, and every accepted state and PDF are saved.

Time halving still changes J by up to 1122.5109 A/m² (at 0.395 us), and E by
660.0122 V/m (at 0.270 us). Population L1 and Te differences are 4.32606e-6 and
0.000259630 eV. The small final current differs by 0.00538928 A/m², about 170%
of the finer tail. Accurate frozen-source oscillations therefore do not establish
resolved reaction/collision splitting or a validated current tail. Further
source/time treatment must retain this limitation rather than infer accuracy
from near-roundoff conservation or small population differences.

The new spatial components and complete longitudinal wave evidence are described
in `ATOMIC_SPATIAL_FIELDS.md`; they retain separate exact-vacuum/interface and
boundary acceptance requirements.
