# Resource guard for the physical species mass audit

Separate immutable dispatcher for `mass` → `species_mass_audit.run`.
It retains the one-second own-process monitor, 8 GiB/256 MiB bounds, final
resource observation and complete optional `--original` replay verification
from `guarded_species_interface`. No original source or resource threshold is
changed. One worker at a time; collect the actual terminal process exit.

The exact-mass vectors remain diagnostics and are never used to evolve a state.
Source seeds explicitly include this guard and its calculation. No unit tests.
