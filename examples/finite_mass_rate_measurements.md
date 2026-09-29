# Complete finite-mass rate measurement protocol

Twenty-one independent cases retain the nonlinear moving reference at n0=1e18
through nr128, original-density controls through nr64, all three compressible
wave branches through nr32, a slow nr64 case and weak-amplitude slow/fast controls.
They measure the complete unsplit RHS with actual mass-flux n_t; manufactured
rates include original species forces/heat while waves are unforced.

Saves initial state, conservative/species/magnetic rates, full E, integrated
momentum/energy fluxes, density rate, velocity representation and independent
exact target arrays. Rejections preserve original inputs and explicit failure.
Source/configuration/parameter/NPZ hashes and iteration residuals are mandatory.
A completed solve is not an accuracy pass; large slow-wave pressure errors
remain in the records. Each case runs separately via --group. No evolution
or formation conclusion follows from this instantaneous protocol alone.
