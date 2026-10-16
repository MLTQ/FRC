# Complete atomic-source thermal evolutions

Run `python -m examples.impact_ionization OUTPUT`. The destination must be new.
This evolves homogeneous fixed-volume direct-ionization source problems, not
the moving coil pulse or a physical gas-discharge preparation calculation.
No unit tests, clipping, species floors or energy restoration are used.

The hydrogen comparison independently eliminates energy using conservation:
Te(x)=[(3/2)*x0*Te0+chi*(x0-x)]/[(3/2)*x]. Then
t(x)=integral_x0^x dx/[N*k(Te(x))*x*(1-x)]. Adaptive scalar quadrature and a
monotone root solve give the reference, independent of the vector RK4 source
integrator. The published rounded13.6eV threshold is used consistently in
this analytic fixture. Three complete20us trajectories measure time convergence.

Three argon cases evolve all19 charge states and electron thermal energy over
3.5us at fixed N=8.961429541213399e19m^-3, matching the initial central nucleus
density of the pulse. This is not the tube-averaged gas fill. Cases begin at
Te=10eV with100% or10% singly-ionized seed; the third begins at10% and receives
20eV per nucleus of explicitly prescribed thermal energy over the first1us.
The heating rate is2*Edep/tp*sin(pi*t/tp)^2 inside that interval, zero outside.
This thermal source does not specify a circuit, coil or preparation mechanism.
400/800/1600 RK4 steps keep identical physical inputs and full duration.

Every state and derived temperature/charge/atomic-energy history is saved.
Number balance and thermal-plus-atomic energy are compared with independently
integrated analytic heating. Report minimum population, fit-range event content,
and event content involving archival threshold discrepancies. No assumption
of physical accuracy follows from numerical convergence or conservation.

Excitation, metastable populations, recombination, radiation, escape, motion,
electron-distribution evolution and momentum transfer remain absent. Results
measure one necessary reaction channel and its energy cost. They must not be
installed into the moving solver as a complete atomic closure. The sourced
data file and its raw inputs are separately identified/copied into the output.
