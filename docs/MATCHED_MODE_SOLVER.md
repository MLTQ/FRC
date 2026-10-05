# Separately refined continuum mode and scalable full equations

This is an extension of STRONG_WALL_DESIGN.md's bounded, periodic, fully ionized
verification problem. It does not supply formation, exterior vacuum, end-loss,
neutral or laboratory wall physics. Hardware remains prescribed waveforms.
Physical masses, actual density, both species and every original wall equation
are unchanged. No density atmosphere, field filtering or energy repair is added.

## Independent accuracy target

The frozen first-radial-Bessel / N=64 axial original-continuum mode has real
growth 8763.749254/s. It supplies only initial data and an independent linear
comparison. The full nonlinear numerical model evolves actual conservative
states at opposite small amplitudes. Their centered physical difference can be
compared with exp(lambda t) times the reference, with independent pressure,
full electric/current, axis/wall and representation errors. The reference is
not a production-space eigenvector and does not force subsequent evolution.

Fluid element count, fluid axial points and both field directions are now
explicitly separate. Rate archives save every actual input/rate coefficient;
large reconstructed physical samples have reproducible dtype/shape/value hashes.
Independent quadrature 16/20 and complete numerical replay remain required.

## Scalable inverse proposals

`ModalWallCoordinates` uses the same fixed mean physical scales but a separable
radial null basis for the toroidal wall trace, with axial points as the inner
index. Normal-wall fluid slots remain inactive. Matching fluid/field axial
counts are needed by the modal inverse, not by the independently refined rates.
Each Fourier mode of the uniform approximate Jacobian has a sparse fluid block
and a small magnetic Schur complement. Every mode, including zero/Nyquist,
is retained. FFT appears only in inverse proposals, never a state filter.

`UniformWallMetric` preconditions induction with the positive physical uniform
current-inertia metric K0 rather than magnetic H alone. For actual n/nref in
[a,b], min(1,1/b) K0 <= K(n) <= max(1,1/a) K0. This changes only the inner
iterative inverse; all variable-density equations and original residual checks
are inherited. The density is never replaced in the physical operator.

Two-stage Gauss still solves Yi=U0+h sum_j Aij F(Yj), with the original full
rounded-stage, endpoint, positivity and induction acceptance. It is fourth
order, without intentional damping of linear imaginary modes. Unresolved phase
and nonlinear energy error must be measured; no conservation repair is applied.

## Nonlinear proposal experiments

All alternatives are opt-in and preserve prior source closures:

- ModalWallGauss: original chord/Newton algorithm with the scalable inverse.
- AcceleratedModalGauss: accepted-iterate multisecant proposals, local to a step.
- ContinuedModalGauss: previous collocation-polynomial guess and retained secants.
- RecycledModalGauss: retains centered Newton Jacobian products instead of tiny
  accepted-iterate differences. Original finite-difference scale and GMRES stay.
- GuardedModalGauss: approximate proposals get one trial and must halve the
  residual or pass. After four iterations it uses secant/Newton only. Original
  Newton backtracking and all physical acceptance bounds are unchanged.

For inverse P and stored pairs S,Y, G=P Y, the proposal is
-P r -(S-G) least_squares(G,P r). Pair normalization and rank truncation affect
only this approximate inverse. They never remove physical modes. Continued
memory is used only for an exact match to the last accepted physical input and
time step, and commits only after a successful endpoint. Rejected steps leave
inputs and previous history unchanged. Slower pilots remain archived; solver
work counts and measured accuracy are separate evidence.

The complete trajectory protocol saves every accepted state/stage and failure,
with full source/runtime identity. The independent audit integrates original
continuous kinetic, separate internal, physical current and magnetic energy,
including harmonic bias cross work, plus particles and entropy. These are
separate from grid invariants and do not close local interface budgets.
