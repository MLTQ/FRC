# Unsplit moving finite-mass derivative

One stage computes conservative bulk mass flux first. That exact flux gives
primal rho_t, mapped dual n_t and dual particle transport; no second continuity
equation is evolved. All mass divisions use mi+me. CurrentThermal supplies
instantaneous geometry and thermal/current coefficients only; its fixed-n
time step is not inserted into a moving split.

MovingCurrent includes actual n_t, regular bulk velocity and curl v. Its full
E enters Poynting work. Shared finite-volume momentum/energy fluxes include
Maxwell and drift stresses, species enthalpy, bulk/current and mixed current
energy transport, with cylindrical hoop terms. Electron U uses actual bulk
particle flux plus beta-scaled entropy-mean current transport. Ion pressure is
recovered from conserved total E; its adiabatic/local-energy consistency must
converge and is not asserted exact at finite resolution.

FiniteMassSources is explicit manufactured/external force/heat and electric
forcing. Its five fluid rates cannot create mass; ionization is unsupported.
The source callback receives time and the evaluated state. No callback is used
in unforced dynamics. Pressure_rates differentiates the full stated energy
partition, including changing 1/n current energy, without modifying the state.

Returned fluxes are integrated SI rates and expose budgets for independent
audits. Inputs never change. Boundaries remain a filled periodic cylinder with
impermeable radial flow, fixed psi and Btheta(R)=0/free Ez. No vacuum, physical
resistivity/heating, shock entropy guarantee or stable time step follows merely
from this derivative. Numerical bulk diffusion is described in finite_mass_flux.
