# Boundary-conforming bounded plasma candidate

StrongWallFluid declares the outer q=n vr/r slots inactive and exactly zero,
in addition to the original dependent axis velocities. Its pressure reaction
keeps the derivative tangent to this physical impermeable-wall space.
StrongWallPlasma assembles every original bulk/current/species term, then
applies the same zero-work normal reaction to the complete momentum rate.
Density, entropies and tangential momenta are unchanged by that constraint.
The old total-power diagnostic remains valid because wq is exactly zero there;
tests must independently recompute it. Reaction magnitude is reported.

The magnetic/current equations use StrongWallInduction with the declared
variational wall curl and separate raw-curl defect. These are explicit new
discrete boundary equations for the previously intended closed wall. Original
classes and source archives are unchanged. No normal-velocity or magnetic state
is clamped after evolution; stage coordinates must satisfy both constraints.
No exterior, laboratory wall/sheath, end-loss or formation physics is supplied.
