# Reactive species momentum and translational energy

`AtomicMoments` composes the joint atomic network with momentum and total
translational energy for every retained heavy level and the electron fluid.
Keeping level moments separately lets newly populated levels inherit motion
without assigning an arbitrary velocity to an initially empty charge state.
This is a homogeneous reaction source, not spatial or electromagnetic coupling.

Physical masses are m_a=m_neutral-q_a*m_e; binding mass defects of order
binding_energy/c² are neglected consistently with nonrelativistic mechanics.
The network's former electron thermal coordinate stores lab translational energy;
the temperature passed to rates subtracts actual bulk kinetic energy first.
Appended moments are P/sqrt(m_neutral*eV) and heavy lab energies in eV per
initial nucleus. Number, charge, total mass, vector momentum and
sum(lab energies)+atomic binding+produced photons-integrated heat are independent
evolved invariants. No reservoir is overwritten to enforce them.

For a->b+e, products inherit the donor velocity distribution: P_b gains
(m_b/m_a)P_a per converted fraction and P_e gains (m_e/m_a)P_a. Translational
energy transfers by the same mass ratios; the existing physical ionization
threshold is removed from electron energy. Excitation and spontaneous decay
transfer the whole heavy translational distribution between equal-mass levels.
These ionization moment assumptions follow the small-drift Maxwellian treatment
in [Meier and Shumlak, Eqs.21–34](https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/Meier_PoP2012.pdf).
We retain electron mass and resolved binding energy rather than their reduced
single-fluid/effective-ionization-energy model.

For radiative capture, independently use two-body centre-of-mass kinematics.
With incoming masses m_i,m_e, product M=m_i+m_e, mean velocities v_i,v_e,
ion random energy u_i and fitted mean captured-electron energy c_e, the product
energy per event is

    |m_i*v_i+m_e*v_e|²/(2*M) + (m_i/M)*u_i + (m_e/M)*c_e.

Produced photons receive binding release plus

    (m_i*m_e/M)*|v_i-v_e|²/2 + (m_e/M)*u_i + (m_i/M)*c_e.

Both formulas use joules; code converts to eV. Their sum equals the full
incoming energy plus binding release. The mean photon momentum and photon
recoil energy are neglected at leading nonrelativistic order, not established
zero. Capture selection uses the existing Maxwellian rate derivative moment;
finite drift and heavy thermal speed must be small relative to electron thermal
speed. Diagnostics report both ratios. Photon production is not escape.

The frozen generator transfers moment coordinates as well as populations,
including into initially empty states. Exponential midpoint predicts ne, Te and
electron velocity; it is a second-order nonlinear step. Full linear conservation
does not guarantee positive thermal energies at unresolved timesteps: invalid
stages/endpoints reject without clipping, population floors or energy repair.
An empty species must have exactly zero momentum and energy.

The small moment generator is already assembled as a dense array (103 by103 for
ground-ladder argon including the affine coordinate). Pass that array directly
to SciPy's existing `expm_multiply` instead of constructing a fresh CSR matrix
for each half/full exponential action. The physical entries, supplied trace,
exponential-midpoint stages and SciPy accuracy policy are unchanged; no series
term, population or coupling is deliberately omitted. Dense and sparse matrix
products can differ in roundoff. All ten complete independent DI/RR and heated
argon histories in `outputs/atomic_dense_exponential` succeed, retaining zero
patterns with state differences at most 3.56e-15. DI/RR time-halving errors
decrease approximately fourfold. The full20ns moving control retains 971 steps,
all cell counts and velocity indices, and Erel 2.84209e-11. Corresponding-cell
temperature differences are at most 2.21e-15eV, E/J relative L2 differences
1.06e-14/9.89e-9. Precise-union pointwise maxima include displaced discontinuities;
the comparison preserves those maxima rather than claiming pointwise equality.
Sources and complete evidence are in `reactive_material_argon_profile_10ns`
and `reactive_material_argon_dense_exp_10ns`. Timings are not a controlled
whole-simulation speed benchmark; unresolved physical accuracy is unchanged.

Velocity reconstruction divides the stored dimensionless momentum by population
before multiplying by its physical scale. This avoids an intermediate kilogram
mass times population underflow; no population or momentum is clipped. Kinetic
energy similarly evaluates energy per particle before multiplying by population.

`FixedFrameAtomicMoments` in atomic_frame.py composes this unchanged source with
a fixed Galilean coordinate choice. Its input and observed budgets are lab-frame;
its inherited `fields` and generator operate on relative moment coordinates.
The separate companion documents the transform and retained finite-range limits.

All original missing excited channels, branching, RR-range, opacity and kinetic
limits remain. Coulomb/neutral scattering, charge exchange, spatial transport,
photon momentum/absorption and self-consistent fields are absent from this source.
The strong pulse's large drift can violate this closure; do not install it there
merely because homogeneous balances close. No provider rate arrays are embedded.
