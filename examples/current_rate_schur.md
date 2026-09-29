# Positive magnetic-space static rate solve

Auxiliary verification solver for a purely toroidal field, fixed positive n,
and the same full/isolated convective sources as current_rate_measurements.
It changes no production time evolution. With electric mass M, inertia I,
toroidal Faraday block C and diagonal toroidal magnetic energy H, eliminate
`jdot=-M^-1 C.T H gdot` from `M E=I jdot+f` and `gdot=C E`:

    [1 + sqrt(H) C M^-1 I M^-1 C.T sqrt(H)] y
        = sqrt(H) C M^-1 f,       y=sqrt(H) gdot.

This operator is symmetric positive definite because I is nonnegative and M
positive definite. Conjugate gradients apply the existing accurate electric
mass solves without assembling a large mixed LU. No time step, diffusion,
mass fit or approximation to the physical inertia is introduced. Density may
vary. Nonzero poloidal fields/forcing reject instead of being silently omitted.
The implementation needs only the existing NumPy/SciPy runtime.

After solving, reconstruct electric/current rates and check the original
Ohm, Ampere and Faraday equations plus the eliminated system. Independently
compare with the original mixed sparse direct solve on small meshes and with
Fourier/Chebyshev physical continuum inversion on all measurement meshes.
