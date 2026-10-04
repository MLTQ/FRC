# Independent local current-enthalpy audit

Reconstructs the archived actual n,C,Cdot,J and weak pressure field with
root-product polynomials/FFT. Differentiates the reconstructed physical EOS h
and the current analytically, and measures the original local identity

    epsilon_s,t - J dot E_pressure,s + div(a_s*gamma_s*h_s*J/e) = 0.

The actual reconstructed strong div(J) is included; neither it nor a global
power residual is subtracted. Reports volume/axis/wall/max errors and integrated
terms at quadratures 16/20, with independent source/raw hashes and a separate
frozen archive. The integrated piecewise divergence does not automatically
include distributional internal-face jumps; these remain part of the later
physical interface-flux gate. This is a smooth local consistency audit, not a
complete moving total-energy or Poynting-boundary proof.
