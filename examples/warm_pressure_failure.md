# High-resolution pressure rejection replay

Reproduces the supplemental nr128 failure at 320/640/1280 steps. Uses exactly
the original physical state, fixed density, pressure midpoint and conservative
face energy flux. It validates the prospective ion partition before assigning
the next accepted arrays. No total energy, pressure or density repair occurs.

Saves lossless initial, last accepted and rejected fields, all accepted histories,
the rejected time/cell/minima/residual, original and actual configurations,
engine and replay/control source hashes. `completed` means the requested final
time was reached; `experiment_complete` means the observation of success/failure
finished. These flags are deliberately distinct. Rejected fields are not
presented as accepted physics. The regular pressure measurement raised before
archiving at nr128; this dedicated replay preserves that failure evidence.
