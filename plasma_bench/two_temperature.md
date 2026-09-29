# Equal-gamma species thermodynamics

`SpeciesThermodynamics` places electron U and pe on the scalar dual grid,
with n equal to the **overlap of the actual ion mass divided by mi**.
For singly charged ions, massless electrons and common gamma:
U=pe/(n/nref)^(gamma-1), pe=(n/nref)^(gamma-1) U.
`pressure_rate` includes both Udot and (gamma-1) pe ndot/n; omitting the latter
would lose compression work. The reference nref changes coordinates only.

Electron cell energy is the Q1 average pe/(gamma-1). Ion cell pressure is
the conservative state's total thermal pressure minus that average. The
overlap/average adjoint identity equates global dual and primal electron
energy exactly. No total-energy update or reset is performed here.
For this reduced smooth ideal model the partition converges to both adiabatic
species equations. At finite resolution all residual numerical thermal energy
goes to the recovered ion component; this is numerical error, not a collisional
heating model or verified shock electron-heating prescription.

Every density, U, pe and pi must be positive; exact zero-temperature electrons
use the existing cold solver, not a log(0) convention. No floors. The single
gamma is explicit and shared; different species gammas are unsupported.
Diagnostics separate species energies, minimum temperatures, conserved linear
U integral and physical log entropy divided by kB/(gamma-1), up to a fixed
reference gauge. Finite upwind/Rusanov diffusion need not conserve the latter.

`AdiabaticSpeciesMHD` calls validation for both trial stages before commit.
The standalone fixed-density `ElectronPressure` algorithm remains independent
until pressure/Hall enthalpy and boundary work have been consistently coupled.
