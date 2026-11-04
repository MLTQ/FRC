# hybrid3d_electrons.py

Legacy experimental entropy closure, retained with
`thermal={'closure':'entropy', ...}` for reproduction of invalidated startup
data. The production thermal default is `EnergyElectrons`; its conservative
energy/compression equation avoids the measured pressure-work defect below.

Evolving thermal electrons for full-shot kinetic-ion runs, instead of a fixed
temperature bath. Scalar pressure p=S n^γ, γ=5/3, uses the actual PIC density.
The continuum equations are D_e S/Dt=(γ−1)ηJ²/n^γ and u_e=u_i−J/(ne).

- `pressure`, `temperature`: adiabatic compression plus advected entropy and
  deposited heat; no imposed temperature floor in pressure or added mass.
- `eta`, `edge_eta`: constant η for independent decay comparisons, or NRL
  transverse Spitzer η=1.03e−4 lnΛ max(Te,1 eV)^−3/2 Ω m, lnΛ=10. The 1 eV
  applicability cutoff is explicit. [NRL Plasma Formulary](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Plasma_Formulary_2023.pdf).
- `joule_density`: edge ηJ² transferred to nodes, preserving its volume sum.
- Bulk entropy travels with the same weighted macro-particles that deposit
  ion density. `begin_step`, `prepare_density` and `finish_step` map it with
  CIC; no independent bulk grid advection smooths entropy against an
  unsmoothed PIC density.
- `advance` applies only relative Hall entropy transport, −J/(ne), by
  trilinear semi-Lagrangian remapping per field substep. It accumulates RK3
  Joule energy for deposition at the end of the ion step.
- Bulk PIC transport precedes relative Hall remapping (first-order Lie
  splitting). Transferring the Hall grid change back to carriers uses
  fractional entropy losses and absolute gains. It preserves positivity and
  the requested weighted entropy increment, avoiding heat manufactured by
  clipping negative additive carrier updates. No energy rescaling is used.
- Heat is distributed with the PIC weights and the discrete energy gradient:
  Ue=Σ_p w_p S_p gather[n^(γ−1)]/(γ−1). Each particle's increment is local;
  no global temperature/energy rescaling is applied. `heat_deposition_residual_j`
  measures roundoff/conservation of that transfer.

Hall remapping and the thermal/field splitting are first order and diffusive;
bulk entropy uses PIC transport and the Yee field step is still RK3. The
declared Ohm-law density cutoff also bounds retained heat.
Heat not retained at n<n_floor is measured as `cutoff_removed_heat_j`; it is
not silently deposited into a numerical atmosphere. Exterior vacuum is not
advanced with this electron model. There is no ionization, conduction,
electron-ion temperature relaxation or sheath model yet.

Zero-density cells carry no particle entropy or heat. Once per ion step,
`prepare_density` locates material nodes; empty-node entropy is extrapolated
solely as a transport ghost value. Discarded heat cannot reappear on refill.
`remap_energy_j`, entropy clipping count and absorbed carrier electron energy
are reported separately; boundary work and pressure-force conservation still
require convergence. The remap counter is evidence, not an energy correction.
An optional initial temperature array supports physical entropy profiles.
Fully periodic fixtures wrap all axes; finite device walls use ghost extension.

`DeviceHybrid` removes ηJ from the ion electric acceleration in this
collisional mode: ion-electron friction cancels that electric term. The full
ηJ remains in Faraday's law and its lost magnetic energy heats electrons.
This preserves the continuum resistive energy exchange rather than heating
both ions and electrons from the same work. A complete discrete energy
budget still requires pressure/particle/boundary verification and convergence.
