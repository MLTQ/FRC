# Conservative numerical thermal diffusion

An opt-in replacement for the discontinuous *numerical* pressure contribution
of first-order entropy upwinding. The central species entropy/pressure work
adjoint, full current inertia/Hall terms, moving fluid, physical collision heat,
source work and exact empty vacuum remain. This is not a model of collisional
or nonlocal physical heat conduction, and does not resolve the physical edge.

For fixed geometry let m_i be actual dual-volume particle counts, s_i entropy,
and h_i the particle-weighted internal energy per particle. The EOS gives
h_i=a_i exp(s_i), with positive geometry-dependent a_i. For species s, the
oriented edge particle flux is c_s q, q=-G D.T b/e. Define numerical conductance
g=|c_s q|/2, incidence I, and L=I.T diag(g) I. The added thermal evolution is

    M dh/dt = -L h,
    ds_i/dt = -(L h)_i/(m_i h_i).

Thus d(sum m_i h_i)/dt=0 exactly in the semidiscrete equations. For common
species gamma, total physical entropy production divided by kB is

    sum_edges g*(h_head-h_tail)^2 /
                ((gamma-1)*h_head*h_tail) >= 0.

The generator -M^-1 L has nonnegative off-diagonals and zero row sums, so its
exact exponential preserves positive h and its constant-state limit. We use
the sparse matrix exponential action; invalid numerical results reject without
clipping. Returning to s uses the EOS log(h_new/h_old), not a fitted energy
correction. The measured residual uses the energy represented by that saved s.

For uniform particles on a straight mesh, the coefficient is the usual
first-order upwind thermal diffusion D_num=|u_s|*dx/2. It vanishes with mesh
spacing for smooth fields. For nonuniform density this diffuses actual h rather
than entropy; it changes the spatial numerical scheme. In particular it need
not preserve a local entropy minimum, although total entropy increases. Neither
energy conservation nor total entropy production guarantees spatial accuracy.

Central entropy transport uses the average trace and its original matched
pressure force. It is continuous through q=0. The added numerical thermal flux
also vanishes continuously there and transfers energy only between thermal
volumes; it adds no sign-dependent pressure/current force. This removes that
specific artificial switching mechanism without deleting a physical force.
Other instability, central-transport errors and physical omissions can remain.

`EnergyDiffusionMaterialCoilCurrent.step` uses diffusion half-steps around full
central RK4. Geometry, q and a are frozen only in each diffusion substep;
the complete central stage evolves all states and actual prescribed sources.
Smooth-time order is at most two, potentially lower at nonsmooth fluxes. Each
half-step has no external or vacuum heat flux; no boundary particle escape is
invented. This does not establish the laboratory applicability of a closed
material inventory, especially after actual wall contact.

Saved split diagnostics distinguish entropy production, zero-net thermal-energy
residual and half the sum of absolute nodal energy changes. The last is a
redistribution measure, not generated heat or an integral of every edge flux.
Numerical bulk-viscous ion heating and physical collision heat remain separate.
All cumulative observations are diagnostic only and never modify a state.

The general conservative-central plus dissipative-flux strategy is discussed by
[Ranocha (2018)](https://arxiv.org/abs/1701.02264). The equations above are the
explicit derivation for this material graph, not a claim that a published Euler
flux proves stability of this two-fluid algorithm. Complete coupled runs and
independent transport limits are required before interpreting a new trajectory.
