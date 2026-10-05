# Original parent rates and continuous energy

Reconstructs all fields with the independent root-product polynomial / Fourier
TransportSampler, without production interpolation or differentiation. Scalar
entropy is (S C)/(S n); transverse velocity uses the declared regular quotient
representation. Physical species pressure and its directional time derivative
come from the actual EOS. Errors are measured against independently analytic
contact or linear-wave rates, with global, first/last-element and maximum norms.
The wave uses perturbation norms; contact's exactly zero acceleration is scaled
by velocity/transit time and zero B/J rates are reported as absolute maxima.

Continuous total power integrates bulk kinetic, both internal energies,
mu J²/(2e² n), and full B²/(2mu0), including local bias cross work. This does not
replace the separate grid identity or establish local enthalpy/Poynting face
conservation. The field current coefficients come from production Ampere, then
are reconstructed independently. Quadrature refinement must be measured.
