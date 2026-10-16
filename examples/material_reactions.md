# Complete reactions from different material supports

Builds a physical-moment material state with neutral/ion/electron half-widths
15/5/10mm, densities(.8,.2,.2)×1e20m^-3, temperatures(.1,.1,30)eV, physical
hydrogen masses, distinct vector velocities and an affine velocity variation.
The initially charge-separated state is intentional. Its fixed geometry and
fields are an isolated reaction-source limit, not a viable unforced microsecond
plasma trajectory. Actual transport/forces must surround this substep eventually.

Uses the existing declared constant-coefficient ionization fixture k=1e-14m³/s,
13.6eV binding, with no recombination/radiation. This is a mathematical source
limit, not an atomic-data fit. Reactants are gathered from their separate
material cells. Outside neutral/electron overlap, no ionization is possible;
inside it, ions can be born outside their old5mm support without a seed floor.

An independent pointwise solution applies even where initial charge is nonzero.
Let A=n_neutral(0), C=n_e(0), K=A+C. Then
n_e(t)=K/[1+(K/C-1)exp(-kKt)] for C>0, with conversion n_e(t)-C. Where C=0,
all populations remain unchanged. Per-event donor distribution/mass ratios and
13.6eV threshold independently give every momentum/energy coordinate. Ion and
electron birth velocities coincide, so ionization leaves charge and current
unchanged despite carrier growth; pre-existing current is never reset.

All25/50/100-step histories cover the same2µs reaction duration; two fixed-frame
choices describe identical lab motion, and2/3-point Gauss orders assess gathering
of the piecewise polynomial initial fields. Every accepted extensive state,
actual geometry, integrated budgets, local current/charge changes and new ion
inventory beyond the old support are saved. Empty exterior and unseeded neutral
regions remain exact, without numerical matter or assigned electrons.

`python -m examples.material_reactions FRESH_OUTPUT` with the usual headless
environment creates JSON, raw histories, frozen source and PNG/PDF. No unit tests.
Results remain in ReactionOverlapState: conservative return to moving transport
with variable masses and potentially newly occupied support is still required.
The source operator alone does not complete that coupling or the full FRC goal.

The completed example output additionally retains `gather_and_frame_comparison.json`.
Its independent input integrals use P=sum[mN(vL+vR)/2] and
K=sum[mN(|vL|²+vL·vR+|vR|²)/6], plus1.5eNT and13.6eN_ion. The same saved
histories supply lab/frame and2/3-point integrated comparisons; no states are
changed or replayed for this postprocessing. Equality of Gauss orders here
reflects this constant-rate polynomial case, not general nonlinear quadrature.
