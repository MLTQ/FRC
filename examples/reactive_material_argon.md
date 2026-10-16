# Full-duration charged reacting material attempt

Requests20ns at32/64/128 initial cells per species. Ion/electron initial conditions
retain the original1e14/m³,10mm half-width, Ti/Te=.2/2eV and physical argon masses.
An explicit neutral component at1e14/m³ and.1eV is added; this is a new reacting
case, not an assertion that the previous no-neutral case is physically unchanged.

The actual existing argon ground-ladder ionization and radiative recombination
data drive all retained charge populations, momenta, threshold energy and produced
photons. Their hashes are saved. No provider-private inputs are copied. Original
ground-state Maxwellian, missing excitation/channels and RR-range limitations
remain. Collisions, magnetic/transverse fields, real walls and kinetic validation
are not silently represented by this electrostatic planar experiment.

Source-half / material-full / source-half uses actual changed masses/supports.
Timestep is explicitly half the current acoustic/plasma bound, with unchanged
physical electron mass. The return projection uses8 domain-scale ULPs for redundant
internal cuts only, preserving all species-support endpoints and integrated moments.
No density floor, hidden retry, omitted species or thermal/energy repair.

Declared resource limits are10000 occupied material cells,5000 accepted steps,
and a minimum computational step2e-17s (one billion such steps for20ns).
The timestep limit detects impractical mesh collapse; it does not increase the
accepted step, discard small cells or count a shortened trajectory as success.
The final smallest cell, its actual population/support, and all occupied species
are recorded to distinguish numerical partition cuts from true support ends.
Exceeding one marks the full requested case FAILED; it neither counts a prefix
as success nor changes the requested20ns. Invalid thermal projection or a new
time restriction likewise records failure and the actual last reaction-overlap
state. Every accepted ragged state, actual time, particle/energy/photon budget,
and source archive is retained. This deliberately challenges the current coupling
beyond its soluble neutral limit; completion/convergence must be measured.

Usage: `python -m examples.reactive_material_argon FRESH_OUTPUT` in the usual
headless environment. No unit tests or generic pass/fail framework is invoked.
