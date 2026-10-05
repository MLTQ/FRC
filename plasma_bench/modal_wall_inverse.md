# Fourier-block implicit proposal

Factors I-hL independently for every periodic axial mode, using sparse fluid
LU and a small radial magnetic Schur complement. Applies forward/inverse FFT
only to correction vectors. Keeps the original complex Gauss stage transform
and supplies real stage corrections. Matching axial counts are enforced by the
coordinate class. Shift changes refactor all blocks; an unchanged shift reuses
factors. Positive finite real shift and explicit factor-before-solve contracts
match the earlier dense inverse.

There is no mode truncation, spectral filter, physical diffusion or residual
relaxation. Construction/storage scale linearly in axial count for fixed radial
resolution. All modes, including zero and Nyquist, participate. Independent full
operator and real two-stage substitution controls are required before use.
