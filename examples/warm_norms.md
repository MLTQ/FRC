# Physical species errors

Independent Gauss integration in x=r²/2,z compares actual continuous fields
with the exact contact/simple-wave reference. Electron pe and mapped n are
Q1 interpolants of their dual coefficients. Total fluid pressure/density use
the production limited reconstruction, explicitly sampled at the quadrature
points. Reconstructed pi=ptotal-pe; Te=pe/(ne e), Ti=pi/(rho/mi e).
This ion reconstruction is a diagnostic convention for the smooth fixtures,
distinct from the stored cell ion pressure. It is not a species trace limiter.

All errors use full cylinder volume. Relative norms divide by exact deviations
from the uniform reference background, not the much larger background itself.
RMS SI errors and minimum reconstructed values are also retained. A zero
reference perturbation returns null. Orders 8 and 12 estimate norm quadrature
uncertainty, separately from discretization and temporal errors.
