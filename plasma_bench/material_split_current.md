# Integrate the entropy/current upstream switch

Purpose: an opt-in time method for the actual moving finite-coil simulation.
The existing upwind entropy flux is continuous at zero current; its energy-paired
pressure force can jump there. Explicit RK stages can repeatedly cross this
switch. This module separates that numerical diffusion from the central physical
entropy/pressure pair and integrates its switch implicitly. Two complete weak-pulse
trajectories now measure its behavior; broader physical cases remain opt-in.

Let S select material dual volumes, m their fixed particle counts, I the oriented
edge incidence, G=diag(edge_scale), q=-G D.T b/e, and c_s the species drift fractions.
The upstream value is s_avg - sign(c_s q)*delta(s)/2. Subtracting central transport
leaves entropy diffusion with graph Laplacian L_s=I.T diag(|c_s q|/2) I. Its
matched magnetic force is D G k sign(q), with

    k = -sum_s |c_s| delta(s_s) delta(h_s)/(2e),
    A = H + mu/e² D Mn^-1 D.T,  K = G D.T A^-1 D G.

`CentralEntropy` uses the average trace for both flux and pressure work; potentials,
heat, masses and every physical term remain those of the original model.
`SplitMaterialCoilCurrent` applies a diffusion half-step, the original RK4 with
that central pair, and another diffusion half-step. Geometry/momenta are fixed
only within the diffusion substeps. Whole-method temporal order is at most two
for smooth solutions, and can reduce at a switching event; no RK4 order claim.

For duration tau, solve q_mid=q_old-tau*K*eta/(2e), eta=k*sigma. The trace sigma
is sign(q_mid) off zero and lies in [-1,1] at zero. This is the bounded completion
of the discontinuous upstream pressure, not a resistivity or entropy floor.
Positive k gives a convex box problem in eta, with a retained working set in
`material_current_box` after the bounded minimization. Newly hit bounds are kept
until the free equations settle, then a wrong-sign constraint can be released. Compatible null gradients are retained via a
least-squares free solve, not regularized with an artificial penalty. Iterations
reuse the previous trace to avoid introducing a new sign flip in the solve.
Negative k retains its incoming
branch and must settle to a consistent midpoint sign. Nonconvergence rejects;
uniqueness of a nonconvex branch is not claimed. All negative terms are retained.

The entropy solve is (M+tau L/2) delta(s)=-tau L*s_old. Its energy gradient is
h_discrete=h_old*expm1(delta(s))/delta(s), with limit h_old. Using s_mid and this
h_discrete in k pairs the exact thermal increment with the quadratic field/current
energy increment. Iteration solves these equations, never a fitted energy budget.
The measured energy exchange residual is reported and never fed back into state.
The central split alone has no entropy minimum principle. Midpoint diffusion also
requires small enough graph outflow for monotonicity; whole-run extrema and
current/time differences must be measured, without clipping any state.

Source response power based on an arbitrarily chosen instantaneous trace at zero
is not an independent integral for this method. Use directly observed -Idot*lambda
for external work, retaining physical induced/current energy and actual Joule heat.
All source/geometry/momentum/transport closures are otherwise inherited. Exact
empty vacuum is unchanged. No new physical edge, ionization or stability claim.

Contracts: optional subclass used by `examples.material_coil_pulse`; saved states
retain MaterialCurrentState format and can be sampled by existing visualization
and field comparison tools. Central-species selection is restored after each
stage, including failure; a model instance is not a concurrent integration object.

Measured full1000/2000-step weak pulses reduce final time-halving J difference
to.38724%; weak Jr/Jz differences3.725%/4.055% replace the original15.28%/14.99%.
This is improved time behavior, not spatial convergence or a strong-pulse result.
Both entropy minima remain within2.92e-13 of initial without clipping. See
`docs/MATERIAL_EXPERIMENT.md` for complete measured results and preserved failures.
