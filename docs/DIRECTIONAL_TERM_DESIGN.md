# Separate resolution, physical cancellation and arithmetic-path diagnosis

## Scope

Continue the fixed radial16,N128 original-continuum study without changing any
accepted solver source. Same filled periodic physical-mass fixture and wall
space as RADIAL_MODE_DESIGN.md. Change axial points, radial field cells and
fluid radial elements independently. Preserve complete signed inputs, accepted
rates or explicit rejection, original strict bounds, and full bitwise replays.
No transport coefficient, numerical diffusion, atmosphere, filter or energy
repair is introduced. All findings are verification evidence, not lab settings.

## Original parent terms

The independent radial_parent_terms reference evaluates full ideal, Hall,
pressure and physical current-inertia E, and their negative cylindrical curls.
For F=(p(z)J1(qr),t(z)J1(qr),w(z)J0(qr)),
-curl F=(t'J1,-(p'+q w)J1,-q t J0). Pressure retains the perturbation of 1/n
multiplying the background pressure gradient. Inertia is mu lambda J/(e² n),
using physical masses and the independently solved growth rate. Pressure and
Lorentz acceleration are independently differentiated as original forces /Mn.
Off-grid fourth-order finite differences check each curl; sums check the
independent original mode E and its complete momentum and induction rates.

Captured production terms use actual signed conservative inputs and original
operators. Bulk + electromagnetic + current stress + species + wall reaction
must reconstruct the archived full fluid rate bitwise. Original full Ohm is
rechecked against saved E and Jt, with exact input hashes. Full and reduced
pressure plus its analytic potential remain separately visible. No term is
zeroed, fitted or removed from an accepted derivative.

## Independent physical audit

Root-polynomial/Fourier reconstruction at quadrature16/20 is batched only in
radial rows to bound memory. Every actual sampled term is streamed into a
fingerprint and per-radial-field-cell actual/reference/error squared integrals.
Report physical component L2 norms, errors normalized by the corresponding
small net reference rate, and own-term relative errors when its reference is
nonzero. Axis/wall norms use the local reference; squared-error fractions
separately establish how much of the global error is carried by those regions.

All pairwise error cross integrals are retained. A sum of error magnitudes
cannot diagnose cancellation between large terms. Current-inertia error uses
the saved numerical Jt and therefore includes the result of the induction
solve; it is not an independent imposed forcing error. Independent strong
volume curl and the explicit wall-projected Faraday representation are distinct.
Strong curls exclude surface distributions; existing one-sided face audits
remain necessary for complete local budgets.

## Failed solves and driving diagnostics

A rejected induction case still has well-defined original forcing and momentum
terms at its saved positive physical state. A separate driving capture/audit
can evaluate these without solving induction. It never supplies or accepts a
magnetic/current derivative and does not invent zero inertia. Every record
retains the original full-case status and full_plasma_rate_accepted=false.
Its driving-sum target is ideal+Hall+pressure, while amplification by the full
continuum B-rate norm is explicitly labeled. This is not a passing B-rate error.
The accepted baseline must match shared coefficients bitwise before using this
control to interpret rejected radial refinement.

## Arithmetic-path observation

The instrumented balanced core leaves original iteration and acceptance logic
unchanged and saves each proposed coordinate vector. Let u=A(y)-rhs and
v=(HT)^T(Fw E-Ty); in exact arithmetic u+v=0. A measured gap larger than twice
the unchanged tolerance means these two residual norms cannot both pass at that
iterate. This alone does not prove every mathematically possible iterate fails.

Measure d=(HT)y-H(Ty), and predict the gap as
T^T(HT)y-(HT)^T Ty + (HT)^T Fw M^-1 I M^-1 Fw^T d.
The remainder measures other rounded mass-inverse/distributive effects; no
coefficient is fitted. The exact saved rejected input must reproduce the
original first-attempt failure message and all residuals. A separate accepted
small case checks bitwise identity of instrumented and original iterations.
This is an observation, not a relaxed residual or adopted solver repair.

The observed matrix-association prediction is too small: about 2.7e-13 versus
5.22e-10. It does not explain the failure. A second independent reconstruction
uses S=M^-1, b=Sf, jp=S Fw^T(HT)y and jq=S Fw^T H(Ty), with wp=I jp and wq=I jq.
The opposing residual gap is decomposed as

    (HT)^T Fw [ S(f-wq)-Sf+S(wq) ]
  + (HT)^T Fw [ S(wp)-S(wq) ]
  + T^T(HT)y-(HT)^T Ty.

Each saved vector is evaluated without rerunning CG. The first term, failure
of floating-point mass-inverse distributivity, has norm5.20e-10 versus an
observed5.22e-10. Its relative E difference is9.50e-16. The complete predicted
vector leaves7.06% of the observed gap unexplained; no claim of exact causal
closure is made. At the first original iterate, split E=Sf-S(wq) passes original
dual/physical Faraday at5.28e-11/1.63e-11 and original Ohm at5.23e-16.

## Opt-in linear-inverse arithmetic

SplitInverseWallCore caches Sf for both original RHS and physical E, then
evaluates E=Sf+S(I Jt). This is algebraically the identical linear equation to
E=S(f+I Jt). Original A, actual variable-density inertia, finite physical masses,
boundary space, all residual definitions/limits, CG budget, three inner and
five outer refinement attempts remain unchanged. No error metric is replaced.
The independent cylindrical polynomial controls check full B/J/E and original
rate work on nr6/12/24; iteration exhaustion must leave inputs unchanged.
Separate full exact-captured-input and signed physical-rate studies with fresh
bitwise replay are required. The frozen full trajectory driver does not adopt
this opt-in arithmetic automatically. Spatial accuracy is a separate gate.

## Next acceptance

Separate forcing errors from conditioning and from failure to solve the
original equations. Continue full spatial, temporal, amplitude, residual,
local-work and conservation verification before adoption into growing-response
trajectories, moving exact vacuum, physical transport/ionization, coil/end
boundaries and prescribed-waveform formation. Axisymmetric topology never
establishes three-dimensional stability or particle confinement.
