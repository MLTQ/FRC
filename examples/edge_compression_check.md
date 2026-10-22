# edge_compression_check.py

An analytic plasma–vacuum edge check for bead frc-ibl.11.9: adiabatic θ-pinch
compression of a magnetized column in the 12" device. The equations,
assumptions and flux-conserving end conditions are in the module docstring.

## Case
- **Column:** radius a0 = 0.12 m, with exact vacuum out to the 0.1524 m wall.
  Hydrogen at n0 = 5·10¹⁹ m⁻³, Ti = Te = 20 eV, internal field
  B_in0 = 0.0412 T.
- **External field:** a uniform solenoid field rises 0.05 → 0.15 T on a
  4 µs cosine ramp, then holds to 6 µs.
- **Numerics:** 20 × 20 × 40 mm cells, 8 ions per cell, 5% floor, energy
  electrons, dt = 2 ns.
- **Analytic final state (quasi-static):** radius ratio a/a0 = 0.571, core field
  ratio 3.064.
- **Ideal energies:** for a lossless, infinitely long column, ion kinetic energy
  rises from 23.7 to 56 J and electron energy from 23.7 to 50 J.
- **Where the analytic model is not exact:** the absorbing tube ends, the
  kinetic ion edge, the breathing oscillation and the ion-rms radius (which a
  halo of spilled ions biases outward). It is accurate to a few percent in
  geometry and is not an exact energy reference.

## Results (outputs/edge_compression_*; held phase 4–6 µs averaged)
| model | radius ratio | core B ratio | source work (J) | numerical removal (J) | of which empty-node / cutoff / projection | electron energy (J) | ion kinetic (J) | load residual (J) |
|---|---|---|---|---|---|---|---|---|
| analytic | 0.571 | 3.064 | – | 0 | – | 50 | 56 | – |
| floor | 0.600 | 3.10 | 225.6 | **181.8** | 105.4 / 76.3 / – | 14.5 | 51.8 | −0.38 |
| mask | 0.593 | 3.21 | 54.1 | 30.9 | 27.0 / 0.1 / 3.8 | 17.3 | 37.7 | +0.02 |
| edge, field only | 0.595 | 3.21 | 51.5 | 28.4 | 25.0 / 0.1 / 3.3 | 17.4 | 37.5 | +0.03 |
| edge + electron edge | 0.590 | 3.21 | 51.7 | 19.5 | 16.2 / 0.0 / 3.2 | 21.3 | 37.5 | −0.62 |

Wall time: 14, 18, 39 and 41 minutes respectively. With the interface floor
capped at 2, the edge model needs about twice the field substeps.

## Reading
- **The geometry passes in every model.** The core field compresses to within
  1–5% of the analytic ratio. The radius agrees to 4–5%, mostly from the ion
  halo.
- **The floor medium fails the energy test.**
  - The empty gas annulus is a 1 eV Spitzer resistor, and the external field
    has to diffuse through it.
  - The source therefore does 4× the work, and 81% of it is deleted: 105 J of
    Joule heat removed on empty nodes, plus 76 J cut off.
  - This matters for any configuration that starts with a vacuum region inside
    the gas, such as a partial fill or a localized puff.
- **A sharp edge removes that.** The mask and edge models need only the work of
  the compression itself.
- **They still lose 38–57% of the work, mainly as electron energy removed at the
  plasma edge.** The electrons keep 17 J instead of reaching 50 J. Cause:
  first-order upwind transport of the tangential diamagnetic drift on the
  staircase edge carries energy into empty nodes, where it is deleted.
  - The electron edge closure blocks that flux and uses the interface density.
    It cuts the removal to 16 J and keeps 21 J of electron energy.
  - The cost is the closure's explicit pressure time-centering error, which
    grows from about 0.03 to 0.6 J.
- **The plasma-side interface floor makes no measurable difference** (28.4
  against 30.9 J). It is limited by CIC noise: 1,600–1,700 plasma nodes still
  sit below the floor.
