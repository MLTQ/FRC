# Positive reconstruction without cell-state repair

The initial second-order method reconstructs total pressure and magnetic field
independently. Even when cell pressure is positive, a face can have
`p_T = pt_reconstructed - B_reconstructed²/(2 mu0) < 0`. That happened for the
independent closed-flux fixture at p0=0.1 Pa on 12/24-cell radial meshes.

## Local method

1. Contract density slopes about the original positive cell density when the
   linear-in-x,z reconstruction would produce a nonpositive face or corner.
   Here x=r²/2. The radial density moment and cell state remain unchanged.
2. Independently reconstruct thermal pressure as a linear-in-x,z polynomial.
   Contract its slopes about the positive cell pressure to obtain positive
   auxiliary traces `p_A` at all face/corner evaluation points.
3. Use one weight `w` per cell in `p = w p_T + (1-w) p_A`. For each evaluation
   point where `p_T < b`, constrain `w <= (p_A-b)/(p_A-p_T)`. The smallest
   such bound across the cell's radial faces, axial faces and four corners
   guarantees admissible evaluation traces. Otherwise w=1 and the original
   method is unchanged.

The bound is `b=1e-10*cell_pressure`, with floating-point slack, not a physical
pressure floor in SI units. It constrains reconstruction only. The small relative
margin prevents roundoff from placing a trace exactly at zero. Density uses the
analogous relative bound. Unit-rescaling tests verify this distinction.

The source term uses the same blended total pressure:
`pt = w pt_T + (1-w)*(p_A + B_reconstructed²/(2 mu0))`.
Its radial integral against dr is analytic for the linear pressure parts;
four-point Gaussian quadrature integrates B² exactly for the implemented
quadratic/cubic radial field polynomials. An independent 20-point quadrature
checks the entire modified source.

This is a reconstruction change, not an added heat source. Conserved cell mass,
momentum, energy and corner flux are untouched. Numerical face fluxes remain
shared between neighboring cells; total mass and energy retain their closed-
domain budgets. CT field evolution retains telescoping face divergence. The
limited local pressures can change the truncation error and artificial motion,
so limiter activation counts and refinement remain part of the evidence.

## Scope of the guarantee

Positive face/corner states are guaranteed for valid cell data within floating-
point range. No theorem is asserted here for arbitrary full multidimensional
MHD updates. Magnetic energy and divergence constraints make that a stronger
condition than positive primitive reconstruction. The solver still validates
both RK stages and rejects invalid states without committing repairs.

The local method is independently derived. It is **not** claimed to implement
the full positive CT schemes analyzed by
[Pang and Wu](https://arxiv.org/abs/2410.05173), or the conservative MHD positivity
analysis by [Wu](https://doi.org/10.1137/18M1168017). Those results motivate keeping
the reconstruction guarantee distinct from a full-update guarantee.

`positivity_reconstruction=False` explicitly retains the unbounded reference.
Diagnostics count density-slope contractions, auxiliary thermal-slope contractions
and actual pressure blends over accepted cell-stages; these are not unique-cell
counts. Auxiliary slope limiting has no effect on the applied pressure when
w=1. Minimum blend weight and total cell-stage count make activity measurable.

## Evidence and next step

Tests verify positive multidimensional traces, unchanged input state, consistent
source integration, bit-for-bit unchanged resolved wave evolution, unit scaling,
conservation/divergence and decreasing limiting with equilibrium refinement.
[Measured equilibrium results](EQUILIBRIUM_VERIFICATION.md) retain enabled and
explicitly disabled runs. [Low-pressure/rarefaction results](POSITIVITY_VERIFICATION.md)
separately measure robustness and accuracy; finishing a case is not proof of
adequate resolution or physical model validity.

Dynamic true-vacuum/coil coupling can now be developed against this bounded
reference core. Its boundary/source-work tests, transport/Hall applicability,
formation convergence and stability/confinement questions remain open.
