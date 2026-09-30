# Independent radial force-term diagnosis

Splits the frozen candidate's initial radial acceleration into thermal pressure
source, kinetic/centrifugal source and conservative convection including the
actual density derivative. Independent physical reconstruction and complex-step
analytic spatial derivatives compare each part with original cylindrical Euler.
All four meshes keep individual force arrays and full inputs. This distinguishes
large mutually canceling errors from genuine component consistency; it does not
change the solver or fit a force to the reference.
