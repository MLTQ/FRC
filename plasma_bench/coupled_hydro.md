# Gas-only conservative subsystem

`gas_from_primitive` and `gas_primitive` explicitly use gas energy (internal plus
kinetic), without magnetic energy. The coupled magnetic subsystem owns induced
field energy over the whole plasma/vacuum domain; adding cell magnetic energy
here would double count it.

`hydro_step` advances annular finite-volume Euler equations with existing limited
reconstruction, Rusanov fluxes and SSP RK2. Two reflected axial ghost layers
provide closed slip ends. Cropping reconstructed physical faces prevents any
periodic communication between chamber ends. Ghosts are not exterior fluid
degrees of freedom. Radial walls/axis use existing reflection and metric sources.
Mass and gas energy are conserved under these walls. Invalid stages raise with
no state repair. The coupled solver supplies the combined time-step bound.
