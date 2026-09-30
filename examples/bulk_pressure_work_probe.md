# Conservative pressure force and kinetic-work diagnostic

Constructs a candidate finite-volume pressure stress from actual reconstructed
density and both continuous nodal entropy fields, with physical face quadrature
and the cylindrical pressure hoop source. It is conservative in axial momentum.
The true bulk kinetic-energy pressure derivative uses the reconstructed
conservative momenta, velocity and momentum-force pullback, not cell-center
velocity. Thermal work comes independently from the verified common continuity
and log-entropy component.

Independent analytic cell-average pressure gradients measure force accuracy.
Whole/axis absolute force norms and generalized energy-work discrepancy remain
separate. A uniform-pressure analytic control exposes initial representation
error rather than resetting pressure; opposite flow exposes work sign. The
pressure candidate is diagnostic only and is NOT a production momentum update,
an energy repair or a proof of a complete conservative moving method.

The helper/companion has a separate archive and hash alongside the frozen bulk
calculation identity. Raw force/face/work arrays are preserved for review.
