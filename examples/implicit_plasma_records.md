# Explicit implicit-study source closure

Recursively follows local Python imports from caller-provided calculation/test
seeds and package initialization. Hashes complete source/runtime and reuses the
established immutable source/companion archive writer. Explicit seeds keep pilot
and final protocol identities distinct without implicitly discovering files
created later. No previous frozen provenance implementation is changed.

For a closure containing precise_geometry, runtime identity also records the
loaded geometry backend and precision; MPFR adds gmpy2, MPFR and GMP versions.
This distinguishes identical Python sources run under different arithmetic.
The backend is captured from the actual imported module, not guessed from an
environment variable after the calculation.
