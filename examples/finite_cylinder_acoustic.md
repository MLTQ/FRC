# Full finite-tube acoustic evolution

Runs one complete standing acoustic period in a sealed 1 m by 15 cm diameter
cylinder, using the actual nonlinear two-entropy fluid and explicit end caps.
Ar+ physical mass, n=1e20/m³, Te/Ti=5/1 eV and gamma=5/3 are comparison inputs,
not recommended lab conditions. Initial density is n0[1+epsilon cos(pi s/L)],
velocity zero and each specific entropy constant. Here s=z+L/2.

The independent linear limit is dn/n0=epsilon cos(k s)cos(k cs t),
vz=epsilon cs sin(k s)sin(k cs t), dp/p0=gamma dn/n0, cs²=gamma p0/(M n0),
k=pi/L. End velocity is zero. The pressure difference exerts a physical force
A[p(left)-p(right)] and total axial momentum is 2 M n0 A epsilon cs sin(k cs t)/k.
This tests finite end response, not periodic wrapping. Nonlinear departures are
expected at finite epsilon; amplitude variation distinguishes them from error.

An independent root-product polynomial reconstructs denser axial observations.
Native and reconstructed energy are recorded separately relative to acoustic
perturbation energy, along with mass, momentum, pressure force and constraint
reaction. Subtraction at very small amplitude can expose floating-point thermal
background cancellation and is not a claim of perturbation-energy precision.
All accepted state arrays, source closure and sampled physical histories remain.
The unchanged RK4 method advances all variables with no splitting or repairs.
No unit tests or new reporting/replay framework. This is only the finite fluid
component; current, magnetic vacuum and finite coils must still be connected.

An explicit mode option sets k=m*pi/L and duration=2L/(m cs). The odd-mode
momentum expression generalizes with factor (1-(-1)^m). Higher modes resolve
spatial end response before the linear comparison reaches its finite-amplitude
error. The initial m=1 run remains frozen separately.

Energy changes also use an algebraically equivalent log1p/expm1 difference of
the physical equation of state, retaining changes in density and each entropy.
This avoids subtracting two large thermal totals. It does not alter a state or
repair energy. Both original absolute totals are still recorded.
