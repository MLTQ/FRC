# iec_recover.py

One concrete recovery of the original incomplete HOPS archives. Does not alter
old states or synthesize particle dynamics. Require unchanged physics-module
hashes and NumPy version, and bitwise-identical initial particle draws.

The source RNG is used only for loading and `WallEmitter.inject`; display RNG
is separate. With constant macro weight, geometry and prescribed wall current,
injection draws do not depend on absorption or fields. Replay the real injection
method into a discarded scratch population, avoiding all particle pushes/Poisson
solves. At every original saved time, require exact cumulative birth counts,
emitted charge and injected energy. Any mismatch fails recovery.

Create complete checkpoints for the last two original frames. Restore the actual
particles and Poisson solution, initial energy and all counters, exact RNG and
controller state. Require exact charge/potential/field-energy reconstruction.
A separate actual saved-interval replay must still match the original endpoint
bitwise before production continuation. The original files remain untouched.
